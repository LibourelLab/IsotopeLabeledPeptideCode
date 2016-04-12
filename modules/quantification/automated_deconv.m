function [amino_acid_mid,fval] = automated_deconv(peptides,iterations,start)
%% automated_deconv takes a set of peptides with unique relative isotope abundances and solves for the underlying amino acid MIDs that best describe the collection peptide measurements
% inputs: 
%        peptides: structure array containing amino acid sequence, charge,
%                  a relative isotope distribution measurement (measured) and an isotope
%                  index for the peptide.  there is an optional weight that allows one to
%                  weight the optimization in favor of fitting highly abundant peptides
%       iterations: # of iterations for the deconvolution
%       start: starting point for the amino acid MIDs (127x1 vector) 
%


%define natural abundance MIDS
na_mid = get_elemental_isotope_abundances();

%% find proper sequence field
seq_id = {'seq', 'Sequence','sequence'};
m = isfield(peptides,seq_id);
seq_id = seq_id(m); seq_id = seq_id{1};


% define the amino aid character id
AA = getAminoAcidLabels();
% get the atomic composition for each amino acid
aa_c = cellfun(@atomiccomp,AA,'uni',0);
% find the number of carbons per amino acid
aa_c = cellfun(@(x) x.C,aa_c,'uni',1);

% set optimizer options (turning off display will speed up code)
options = optimoptions('fmincon','Display','iter','TolFun',1e-11,'TolX',1e-11, ...
                    'MaxFunEvals',1e30,'MaxIter',1e4,'Algorithm','interior-point','DiffMinChange',1e-4);

% determine the amino acid distribution per peptide                
[amino_acids]=cellfun(@aacount,{peptides.(seq_id)},'uni',0);


% call equality constraints and set bounds for constrained optimization
A = getLinearEqualityConstraint();
ub = ones(sum(aa_c+1),1);
lb = zeros(sum(aa_c+1),1);
b = ones(length(aa_c),1);

% get the hetero atom MIDS (these are the atoms that do not contain
% variable labeling patterns)
hMID = cellfun(@(x,y) get_hetero_MID(x,y), {peptides.(seq_id)}, {peptides.charge},'uni',0);


% for each defined iteration
for m = 1:iterations
    % use the starting point provided by the user or generate a random
    % point
    if isempty(start)
        x0 = getRandomEstimate();
    else
        
        if iscell(start)
        x0 = cell2mat(start)';
        else
        x0 = start;
        end
    end
    % perform gradient-based minimization to find amino acid MIDs that give
    % the best fit to the measured
    [x{m},fval(m),exitflag{m}] = fmincon(@objective_function,x0,[],[],A,b,lb,ub,[],options);
end

output.aaMIDs = cellfun(@(y) mat2cell(y,aa_c+1),x,'uni',0);
output.fval = fval;

% take best estimate from deconvolutions
[~,k] = min(fval);
fval = fval(k);

% convert amino acid vector to cell array
amino_acid_mid = amino_acid_vector_to_cell(x{k});


%% SUBROUTINES %%

% discrete convolution of peptide: amino acid MIDS --> peptide MIDS
    function [MID] = amino_acid_conv(aaMID,pep_aa)
        MID = 1;
        % for each amino acid
        
        for ii =1:length(AA)
            % for the number of amino acid ii in peptide
            for j =1:pep_aa.(AA{ii})
                % convolve amino acid onto growing MID
                MID = conv(MID,aaMID{ii});
            end
        end
        
    end

% objective function for amino acid deconvolution

    function [residuum] = objective_function(x)
        % conver AA MID to an a
        x_aa = mat2cell(x,aa_c+1);
        expected = cellfun(@(x_r) amino_acid_conv(x_aa,x_r),amino_acids,'uni',0);
        % append hetero MIDs
        expected = cellfun(@(x,y) conv(x,y),expected,hMID,'uni',0);
        
        % reduce expected to measured isotope index
        expected = cellfun(@(x,y) x(y+1), expected,{peptides.isotope},'uni',0);
        % renormalize expected
        expected = cellfun(@(x) x./(sum(x)),expected,'uni',0);
        
        % if there is a weight, compute the weighted fractional
        % disagreement (l1-norm) or do not
        if ~isfield(peptides,'weight')
            residuum = sum(cellfun(@(x,y) norm(x-y,1), {peptides.measured}, expected,'uni',1));
        else
            residuum = sum(cellfun(@(x,y,z) z*norm((x-y),1), {peptides.measured}, expected,{peptides.weight},'uni',1));
        end
        
    end

% generate a random amino acid MID that is feasible
    function [random_aa] = getRandomEstimate()
        random_aa = cellfun(@(x) rand(1,x), num2cell(aa_c+1),'uni',0);
        random_aa = cellfun(@(x) x./(sum(x)),random_aa,'uni',0);
        random_aa = cell2mat(random_aa)';
       
    end

% get the linear equality constraints for each amino acid MID (this ensures
% that the amino acid MIDs sum to 1 for each amino acid
    function [A_s] = getLinearEqualityConstraint()
        A_s = zeros(length(aa_c),sum(aa_c+1));
        g = [0,cumsum(aa_c+1)];
        
        for i =1:length(aa_c)
            A_s(i,(g(i)+1):g(i+1)) = ones(1,aa_c(i)+1);
        end
    end

% simulation the pepetide MID from a an amino acid MID
    function [expected_k] = simulate_peptide_MID(x)
        x_aa = mat2cell(x,aa_c+1);
        % convolve amino acid MIDs to get carbon skeleton distribution
        expected_k = cellfun(@(x_r) amino_acid_conv(x_aa,x_r),amino_acids,'uni',0);
        % convolve hetero MIDs
        expected_k = cellfun(@(x,y) conv(x,y),expected_k,hMID,'uni',0);
        
    end

% simulate the hetero MID (this should only have to computed once)
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