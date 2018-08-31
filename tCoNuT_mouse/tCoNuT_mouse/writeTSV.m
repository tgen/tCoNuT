function writeTSV(fileName,cnaMat,col3Header)
% writeTSV(fileName,cnaMat,col3Header)
%
%   Write TSV files for tCoNuT.
%
%   INPUT:
%       fileName is the name of the file to be written. Desired extension
%       needs to be included. Typical choice would be .tsv.
%
%       cnaMat is M x 3 matrix with the following columns: chromosome,
%       position and col3Header
%
%       col3Header is name of the values in the third column of cnaMat.
%       This name will be used as the column header.  Examples are
%       Fold-Change or BAF.
%

%   [2010] - [2016] Translational Genomics Research Institute (TGen)
%   All Rights Reserved.
%
%   Major Contributor(s):
%       Jessica Aldrich
%   Minor Contributor(s):


fid = fopen(fileName,'w+');

fprintf(fid,'%s\t%s\t%s\n','Chr','Position',col3Header);
for i=1:size(cnaMat,1)
    fprintf(fid,'%d\t%d\t%f\n',cnaMat(i,:));
end

fclose(fid);