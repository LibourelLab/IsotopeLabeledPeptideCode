function [ maxquantProt ] = mapPrecursorMZtoProteinStruct(maxquantProt,mzxmlfile)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

mzxmldata_unlabeled = mzxmlread(mzxmlfile);

% maxquantProt from get PeptidesFrom MaxQuantOutput
% mzxmldata_unlabeled from mzxmlread
for i =1:length(maxquantProt)
    for j = 1:length(maxquantProt(i).peptide)
        maxquantProt(i).peptide(j).mz_mzxml = mzxmldata_unlabeled.scan(maxquantProt(i).peptide(j).scannum).precursorMz.value;
    end
end


end

