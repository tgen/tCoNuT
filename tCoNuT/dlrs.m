function dlrsStat = dlrs(log2Ratios)
% dlrsStat = dlrs(log2Ratios)
%
%   Derivative log-ratio spread is a measure of the noise in the log2
%   ratio data.  
%
%   INPUT:
%       log2Ratios is a vector of log2ratios.
%
%   OUTPUT:
%       dlrsStat is a vector DLRS values.
%

%   [2010] - [2016] Translational Genomics Research Institute (TGen)
%   All Rights Reserved.
%
%   Major Contributor(s):
%       Jessica Aldrich
%   Minor Contributor(s):

xdiffSTD = std(diff(log2Ratios));
dlrsStat = xdiffSTD/sqrt(2);
