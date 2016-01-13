function [ protein ] = getTotalIonsForPro(protein)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

for i =1:length(protein)
    for j = 1:length(protein(i).peptide)
        protein(i).peptide(j).total_ions = sum([protein(i).peptide(j).envelope.ions(:)]);
    end
end

end

