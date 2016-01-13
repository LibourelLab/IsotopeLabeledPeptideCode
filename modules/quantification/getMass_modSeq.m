function [ molweight ,atoms,atoms_no_mods] = getMass_modSeq(peptide)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if isfield(peptide,'seq')
atoms = atomiccomp(peptide.seq);
elseif isfield(peptide,'Sequence')
atoms = atomiccomp(peptide.Sequence);
end

if isfield(peptide,'modseq')
modified_sequence = peptide.modseq;
elseif isfield(peptide,'ModifiedSequence')

modified_sequence = peptide.ModifiedSequence;

end


if isfield(peptide,'charge')
charge = peptide.charge;
elseif isfield(peptide,'Charge')

charge = peptide.Charge;

end


atoms_no_mods = atoms;


expression = '\(\w*\)';
mods=regexp(modified_sequence,expression,'match');
for i =1:length(mods)
    switch mods{i}
        case '(de)'
            % loss of one amine(NH3) and loss of a water moleule, resulting
            % in the net loss of one N and H, and net gain of one O
            atoms.N = atoms.N - 1;
            atoms.H = atoms.H - 1;
            atoms.O = atoms.O + 1;
        case '(ca)'
            % loss of one hyrogen with teh addition of pme
            % carboxymethylated group with iodoacetamide: C2H4INO (loss of
            % one H - I + C2H4INO
            atoms.C = atoms.C +2;
            atoms.H = atoms.H +3;
            atoms.O = atoms.O +1;
            atoms.N = atoms.N +1;
        case '(ox)'
            % assume only one oxidation can occur on methionin
            % this will 
            % need to be changed if we allow double oxidation of methionine
            % to methionine sulfone
            atoms.O = atoms.O +1;
    end
end


    atoms.H = atoms.H + charge;
    atoms_no_mods.H = atoms_no_mods.H + charge;
    
    load constants
    molweight = atoms.C*element.carbon.mass + atoms.O*element.oxygen.mass + atoms.S*element.sulfur.mass+ atoms.N*element.nitrogen.mass + atoms.H*element.hydrogen.mass;

end

