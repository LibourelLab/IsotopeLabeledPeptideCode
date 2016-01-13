function [x_aa] = amino_acid_vector_to_cell(x)

[m,n] = size(x);
if n>m
    x = x';
end


AA = getAminoAcidLabels();
% compute # carbons

% get the atomic composition for each amino acid
aa_c = cellfun(@atomiccomp,AA,'uni',0);
% find the number of carbons per amino acid
aa_c = cellfun(@(x) x.C,aa_c,'uni',1);


x_aa = cellfun(@transpose,mat2cell(x,aa_c+1),'uni',0)';



end