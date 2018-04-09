Data line information for test.VCF:
POS 13273 : Include with CAF above 0.05 and should be printed into baf.txt and merged.vcf.txt
POS 1366274: Include with CAF below 0.05 (Should not be printed)
POS 1339515: Include with CAF below 0.05 (Should not be printed)
POS 10616: No Skip in parseMergeVCF.py because line[10] is a './.' not '.' for tumor and has CAF below 0.05 (Should not be printed)
POS 1361736: No GT or depth information (Should not be printed and throw logging error)
POS 1344411: Test try/catch because of '.' in place of tumor (Should not be printed)
POS 1357461: Include with CAF above 0.05 but sum con_allele and sum aff_allele is below 10
POS 135982: Include with CAF above 0.05 but sum con_allele is below 10
POS 526736: Include with CAF above 0.05 and sum of con_allele and aff_allele is greater than 10 but no rs ID
POS 564598: Include with CAF above 0.05 but sum con_allele and sum aff_allele is below 10
POS 1366433: Normal not called
POS 897730: Should be printed by first if statement into baf.txt and merged.vcf.txt
POS 930567: Should enter 2nd output of if statement and should only be outputted into the merged.vcf.txt Sum of con_allele and aff_allele is greater than 10 but less than 50.
            It also has an rs ID.
