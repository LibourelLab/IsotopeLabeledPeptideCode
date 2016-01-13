function [iscon] = isProteinAContaminant(protein)
% determine if a protein is a contaminant or not

iscon = any(cellfun(@any,regexp(protein.an,'CON\_'))); 

end