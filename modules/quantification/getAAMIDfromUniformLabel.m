function [aa_mids] = getAAMIDfromUniformLabel(label,withHetero)
% define natural abundance MIDS for each element
na_mid = get_elemental_isotope_abundances();



AA = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};


% get the atomic composition for each amino acid
aa_c = cellfun(@atomiccomp,AA,'uni',0);
% find the number of carbons per amino acid
aa_c = cellfun(@(x) x.C,aa_c,'uni',0);
if ~withHetero
    aa_mids = cellfun(@(x) uniform_conv(x,label,[]),aa_c,'uni',0);
else
    hMID = cellfun(@(x) get_hetero_MID(x,0),AA,'uni',0);
    aa_mids = cellfun(@(x,y) uniform_conv(x,label,y),aa_c,hMID,'uni',0);
end


    function [MID,carbonMID] = uniform_conv(C,carbon_label,heteroMID)

        carbonMID = binopdf(0:C,C,carbon_label);
        if ~isempty(heteroMID)
            
            MID = conv(carbonMID,heteroMID);
        else
            MID = carbonMID;
        end
        

    end

    function [heteroMID,c] = get_hetero_MID(sequence,charge)
        atom_types = {'N','O','S','H'};
        atoms = atomiccomp(sequence);
        atoms.H = atoms.H + charge;
        c = atoms.C; 
        
        heteroMID = 1;
        for i =1:length(atom_types)
            numAtoms = atoms.(atom_types{i});
            for k = 1:numAtoms
                heteroMID = conv(heteroMID,na_mid.(atom_types{i}));
            end
        end
    end

end
    