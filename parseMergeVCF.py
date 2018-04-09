# This script parses a SnpEff/SnpSift annotated HaplotypeCaller VCF run on
# control and affected/tumor samples together. It identifies heterzygous SNPs
# in the control sample and grabs the corresponding read counts and allele frequencies
# from both control and affected/tumor samples at the SNP position. In order to filter
# for the high quality common SNPs, the VCF needs to be annotated with dbSNP GMAF/CAF
# which can be done using SnpSift. This script was written specifically for HaplotypeCaller/SnpEff
# but may work on VCFs from other variant callers.

# gatk's HaploType Caller https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
# SnpEff/SnpSift http://snpeff.sourceforge.net/
#
# parseMergeVCF.py norm_gtal-tumor.hc.snpEff.vcf norm_gtalSampleName tumorSampleName readDepth
# Script is python3 compatible

# ARGUMENTS:
#   [0] HaploType Caller SnpEff annotated VCF
#   [1] control sample name used in HaplotypeCaller
#      [2] affected/tumor sample name used in HaplotypeCaller
#   [3] read depth is the min depth of SNP considered for BAF

# OUTPUT:
#   norm_gtal-tumor.baf.txt is text file with header: Chromosome Position NormalReadDepth NormalRefAllele NormalAltAllele TumorReadDepth TumorRefAllele TumorAltAllele CONTROLBAF CONTROLMAF BAF MAF ABSBAFDEV
#      merged.vcf.txt is a text file with header: Chromosome Position NormalReadDepth NormalRefAllele NormalAltAllele TumorReadDepth TumorRefAllele TumorAltAllele

# *  [2010] - [2016] Translational Genomics Research Institute (TGen)
# *  All Rights Reserved.
# *
# * Major Contributor(s):
#    Jessica Aldrich
# * Minor Contributor(s):
import sys
import re
import logging

log = logging.getLogger(__name__)

#Function to assign arguments to the output paths and output files and parse the inputed VCF
def main(vcf_path, normal_name, tumor_name, read_depth):
    out_path = "merged.vcf.txt"
    out_path2 = "{}-{}.baf.txt".format(normal_name, tumor_name)
    string_template = "{split_line_pl1}\t{split_line_pl2}\t{control_pl}\t{con_allele_pl1}\t{con_allele_pl2}\t{affect" \
                      "ed_pl}\t{aff_allele_pl1}\t{aff_allele_pl2}"
    string_template_OFILE2 = "\t{control_baf_pl}\t{control_maf_pl}\t{baf_pl}\t{maf_pl}\t{bafdev_pl}\n"

    with open(vcf_path, "r") as fVCF, open(out_path, "w") as fOFILE, open(
            out_path2, "w") as fOFILE2:  # Double check this with Ryan

        #Write headers for columns in the norm gtal-tumor.baf.txt and merged.vcf.txt
        fOFILE.write("Chromosome\tPosition\tNormalReadDepth\tNormalRefAllele\tNo"
                    "rmalAltAllele\tTumorReadDepth\tTumorRefAllele\tTumorAltAll"
                    "ele\n")

        fOFILE2.write("Chromosome\tPosition\tNormalReadDepth\tNormalRefAllele\tN"
                     "ormalAltAllele\tTumorReadDepth\tTumorRefAllele\tTumorAlt"
                     "Allele\tCONTROLBAF\tCONTROLMAF\tBAF\tMAF\tABSBAFDEV\n")

        #Skip past the meta-information lines of VCF
        for line in fVCF:
            line = line.rstrip('\n')
            if line.startswith('##'):
                continue
            #Retrieve and store index of normal and tumor in header of VCF
            if line.startswith('#CHROM'):
                tmp = line.split("\t")
                norm_gt_index = tmp.index(normal_name)
                tumor_gt_index = tmp.index(tumor_name)
                continue

            split_line = line.split("\t")

            #Check if normal or tumor genotype is not there. Checking using string equality, not containing string
            if ('.' == split_line[9]) or ('.' == split_line[10]):
                log.warning('Skipped line for missing gt: {}'.format(line))
                continue

            #Split Line 8 (GT:AD:DP:GQ:PL) to get index of each variable
            COL8 = split_line[8].split(':')

            #Store index of 'AD' or 'DP' into variables
            try:
                ADind = COL8.index('AD')
                DPind = COL8.index('DP')
            except ValueError:
                continue

            #Split the index line for the normal and tumor genotype and allele depth information
            control = split_line[norm_gt_index].split(':')
            affected = split_line[tumor_gt_index].split(':')

            #Make a list of integers for the alleles for the control and tumor based the index of allele depth
            try:
                con_allele = [int(i) for i in control[ADind].split(',')]
                aff_allele = [int(i) for i in affected[ADind].split(',')]
            except IndexError:
                #Log if there is a missing genotype
                msg = "{} {} {}".format(split_line[norm_gt_index], control, str(ADind))
                logging.warning(msg)
                continue

            # Retrieve the minor allele frequency that is in VCF and store it as the GMAF
            if ('GMAF=' in split_line[7] or 'CAF=' in split_line[7]):
                if('GMAF=' in split_line[7]):
                    GMAF = re.match(r'GMAF=(.*?);', split_line[7]).groups()[0]

                #If its CAF then get the smallest allele frequency
                if ('CAF=' in split_line[7]):
                    result = re.search(r'CAF=(.*?);', split_line[7])
                    if result:
                        CAF = result.groups()[0].split(',')
                    GMAF = min(CAF)

                #Retrieve length of REF and ALT Allele
                charR = len(split_line[3])
                charA = len(split_line[4])

              # Write out baf.txt using a strict filter
                if (float(GMAF) > 0.05
                    #Check if control allele is first allele listed in ALT
                    and control[0] == '0/1'
                    and charR < 2
                    and charA < 2
                    and sum(con_allele) > readDepth
                    and sum(aff_allele) > readDepth
                    and (int(con_allele[1]) / sum(con_allele) > 0.05)):

                    baf = aff_allele[1] / sum(aff_allele)
                    maf = min(aff_allele) / sum(aff_allele)
                    control_baf = con_allele[1] / sum(con_allele)
                    control_maf = min(con_allele) / sum(con_allele)
                    bafdev = abs(0.5 - baf)

                    #Check if chromosome is a sex chromosome
                    if split_line[0] == 'X':
                        split_line[0] = 23
                    if split_line[0] == 'Y':
                        split_line[0] = 24

                    output_format = string_template.format(
                        split_line_pl1=split_line[0],
                        split_line_pl2=split_line[1],
                        control_pl=control[DPind],
                        con_allele_pl1=con_allele[0],
                        con_allele_pl2=con_allele[1],
                        affected_pl=affected[DPind],
                        aff_allele_pl1=aff_allele[0],
                        aff_allele_pl2=aff_allele[1])

                    output_format_f2 = string_template_OFILE2.format(
                        control_baf_pl=control_baf,
                        control_maf_pl=control_maf,
                        baf_pl=baf,
                        maf_pl=maf,
                        bafdev_pl=bafdev)

                    fOFILE2.write(output_format + output_format_f2)

                    # Write out merged.vcf.txt using a loose filter
                if (float(GMAF) > 0.05
                    and split_line[2].startswith('rs')
                    and control[0] == '0/1'
                    and sum(con_allele) > 10
                    and sum(aff_allele) > 10):
                        # Use same output_format as fOFILE2
                        output_format_f1 = string_template.format(
                            split_line_pl1=split_line[0],
                            split_line_pl2=split_line[1],
                            control_pl=control[DPind],
                            con_allele_pl1=con_allele[0],
                            con_allele_pl2=con_allele[1],
                            affected_pl=affected[DPind],
                            aff_allele_pl1=aff_allele[0],
                            aff_allele_pl2=aff_allele[1])
                        fOFILE.write(output_format_f1 + "\n")

if __name__ == "__main__":
    logging.basicConfig()

    vcf = sys.argv[1]
    normalName = sys.argv[2]
    tumorName = sys.argv[3]
    readDepth = int(sys.argv[4])

    main(
        vcf_path=vcf,
        normal_name=normalName,
        tumor_name=tumorName,
        read_depth=readDepth
    )
	

     

