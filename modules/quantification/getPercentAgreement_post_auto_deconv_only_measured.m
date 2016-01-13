function [peptides] = getPercentAgreement_post_auto_deconv_only_measured(peptides,AAMD,flag)
% define natural abundance mids
% na_mid.C = [0.989300000000000,0.0107000000000000];
% na_mid.N = [0.996320000000000,0.00368000000000000];
% na_mid.O = [0.997570000000000,0.000380000000000000,0.00205000000000000];
% na_mid.S = [0.949300000000000,0.00760000000000000,0.0429000000000000,0,0.000200000000000000];
% na_mid.H = [0.999885000000000,0.000115000000000000];

seq_id = {'seq', 'Sequence','sequence'};
fit_id = {'ria','fit','fitting'};
env_id = {'feature','env','envelope'};
m = isfield(peptides,seq_id);
seq_id = seq_id(m);

m = isfield(peptides,fit_id);
fit_id = fit_id(m);
m = isfield(peptides,env_id);
env_id = env_id(m);



AA = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
aa_c = cellfun(@atomiccomp,AA,'uni',0);
aa_c = cellfun(@(x) x.C,aa_c,'uni',1);
options = optimoptions('fmincon','Display','off','TolFun',1e-11,'TolX',1e-11, ...
                    'MaxFunEvals',1e30,'MaxIter',1e4,'Algorithm','interior-point','DiffMinChange',1e-4);

[amino_acids]=cellfun(@aacount,{peptides.(seq_id{1})},'uni',0);

% take best estimate 


exp = simulate_peptide_MID(AAMD);


if ~flag
for i=1:length(peptides)
    peptides(i).simulated = exp{i}';
    sim = peptides(i).simulated;
    meas = peptides(i).measured;
    sim = sim(logical(meas));
    meas = meas(logical(meas));
    
    sim = sim./sum(sim);
    meas = meas./sum(meas);

    
    peptides(i).pa = 1-0.5*norm(sim - meas,1);
   
end

else

[~,c_hMID]=getPeptideMIDFromCarbonLabel(peptides,0.0107);
idx = cellfun(@(x) cumsum(x) > 0.99 ,c_hMID,'uni',0);

for i=1:length(peptides)
    peptides(i).simulated = exp{i}';
    
    sim = peptides(i).simulated(idx{i});
    meas = peptides(i).measured(idx{i});

    sim = sim(logical(meas));
    meas = meas(logical(meas));


    sim = sim./sum(sim);
    meas = meas./sum(meas);

    peptides(i).pa = 1-0.5*norm(meas - sim,1);
   
end
end



    function [MID] = amino_acid_conv(aaMID,pep_aa)
        %amino_acids = aacount(sequence);
        MID = 1;
        for ii =1:length(AA)
            for j =1:pep_aa.(AA{ii})
                MID = conv(MID,aaMID{ii});
            end
        end
        
    end



    function [expected_k] = simulate_peptide_MID(x)
        x_aa = mat2cell(x,aa_c+1);
        expected_k = cellfun(@(x_r) amino_acid_conv(x_aa,x_r),amino_acids,'uni',0);
    end

end