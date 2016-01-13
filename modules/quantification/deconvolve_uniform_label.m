function [percent_label] = deconvolve_uniform_label(peptides,iterations,start)
% fields from peptide are:
% 1) measured ; measured RIA
% 2) isotope: isotope index
% 3) weight (optional)


% define natural abundance MIDS for each element
na_mid = get_elemental_isotope_abundances();

% get heteroMIDs and the number of carbons per peptide
[hMID,numCarbons]= cellfun(@(x,y)  get_hetero_MID(x,y),{peptides.seq},{peptides.charge},'uni',0);

% set optimizer options
options = optimoptions('fmincon','Display','iter','TolFun',1e-11,'TolX',1e-11, ...
                    'MaxFunEvals',1e30,'MaxIter',1e4,'Algorithm','interior-point','DiffMinChange',1e-4);


%define objective function for optimization ------> l1-norm minimization of
%total disagreement between measured and expected RIA

for m = 1:iterations
    if isempty(start)
        x0 = rand(1);
    else
        x0 = start;
    end
    
    [x{m},fval(m),exitflag{m}] = fmincon(@objective_function,x0,[],[],[],[],0,1,[],options);
end


% choose best of all iterations
[~,m] = min(fval);
percent_label = x{m};


    function [residuum] = objective_function(label)
        expected = cellfun(@(numc,hetero) uniform_conv(numc,label,hetero),numCarbons,hMID,'uni',0);
        % reduce expected to measured isotope index
        expected = cellfun(@(x,y) x(y+1), expected,{peptides.isotope},'uni',0);
        % renormalize expected
        expected = cellfun(@(x) x./(sum(x)),expected,'uni',0);
       
        if ~isfield(peptides,'weight')
            residuum = sum(cellfun(@(x,y) norm(x-y,1), {peptides.measured}, expected,'uni',1));
        else
            residuum = sum(cellfun(@(x,y,z) z*norm((x-y),1), {peptides.measured}, expected,{peptides.weight},'uni',1));
        end
        
    end

    function [MID,carbonMID] = uniform_conv(C,carbon_label,heteroMID)

        carbonMID = binopdf(0:C,C,carbon_label);
        MID = conv(carbonMID,heteroMID);

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