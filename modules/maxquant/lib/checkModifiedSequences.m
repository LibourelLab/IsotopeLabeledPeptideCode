function [ peptide ] = checkModifiedSequences(peptide )
%checkModifiedSequences Summary of this function goes here
%   Detailed explanation goes here

mass_diff= peptide.mass - getMass_modSeq(peptide);
modseq = peptide.modseq;
modseq(1) = [];
modseq(end) = [];
%check for extra c-modification
k = round(mass_diff/57);
if k>0
    peptide.mods = strcat('Modified',num2str(k),'(ca)');
    for i =1:k
        modseq = strcat(modseq,'(ca)');
    end
    
end
modseq = strcat('_',modseq,'_');
    
 peptide.modseq = modseq;

end

