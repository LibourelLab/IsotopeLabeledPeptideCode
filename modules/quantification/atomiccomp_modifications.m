function [ atoms,atoms_no_mods] = atomiccomp_modifications(modified_sequence)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

seq = cell2mat(regexp(modified_sequence,'\(\w*\)','split'));
seq([1,end]) = [];

atoms = atomiccomp(seq);

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
            % assume only one oxidation can occur on methionine
            % this will 
            % need to be changed if we allow double oxidation of methionine
            % to methionine sulfone
            atoms.O = atoms.O +1;
    end
end





end

