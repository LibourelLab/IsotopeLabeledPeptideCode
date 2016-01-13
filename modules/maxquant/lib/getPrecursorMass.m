function [ maxquantProt ] = getPrecursorMass(maxquantProt,mzxmlfile)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

if ischar(mzxmlfile)
mzxmldata_unlabeled = mzxmlread(mzxmlfile);
else
    mzxmldata_unlabeled = mzxmlfile;
end
% maxquantProt from get PeptidesFrom MaxQuantOutput
% mzxmldata_unlabeled from mzxmlread
for i =1:length(maxquantProt)
    for j = 1:length(maxquantProt(i).peptide)
        maxquantProt(i).peptide(j).mzRef = mzxmldata_unlabeled.scan(maxquantProt(i).peptide(j).scannum).precursorMz.value;
        maxquantProt(i).peptide(j) = findIsotopeOfReferenceMZ( maxquantProt(i).peptide(j));
    end
end



end

