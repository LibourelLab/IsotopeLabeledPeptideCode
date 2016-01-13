

% convert evidence file into peptide database

%previously identified betaconglycinin peptides
load('beta_con_glycinin_peptides','peptides')
load('/Users/joshuagoldford/Dropbox_personal/Dropbox/software_validation/soy_chymo/peptide identification/maxquant/second run/txt/msms.mat')
load('/Users/joshuagoldford/Dropbox_personal/Dropbox/software_validation/soy_chymo/peptide identification/maxquant/second run/txt/evidence.mat')
load('/Users/joshuagoldford/Dropbox_personal/Dropbox/software_validation/soy_chymo/openms/80813/mass_traces_noPP.mat')

mass_traces = mass_traces_noPP;
clear mass_traces_noPP

e_ids = cellfun(@(x,y) strcat(x(2:end-1),num2str(y)),{evidence.ModifiedSequence},{evidence.Charge},'uni',0)';
pep_ids = cellfun(@(x,y) strcat(x(2:end-1),num2str(y)),{peptides.modseq},{peptides.charge},'uni',0)';
[~,k] = unique(e_ids);
evidence = evidence(k);
e_ids = e_ids(k);

[~,e,l]=intersect(e_ids,pep_ids);
pep_ev = evidence(e);
peptides =peptides(l);

clearvars -except pep_ev peptides mass_traces
for i =1:length(peptides)
    peptides(i).RT = pep_ev(i).RetentionTime;
     u = pep_ev(i).Uncalibrated_CalibratedM_z_ppm_;
     mz = pep_ev(i).m_z;
     k = cellfun(@(x) str2num(x) + 1, strsplit(pep_ev(i).MS_MSIDs,';'));
     %mz = [msms(k).m_z] -
     %[msms(k).IsotopeIndex].*1.003355/peptides(i).charge
     peptides(i).mzRef = mz*(1e-6)*u + mz;
     peptides(i).evidence = pep_ev(i);  

end


bc.peptide = peptides;

clearvars -except bc mass_traces

mtf_51 = filterMassTracesForPro(mass_traces,bc,5,1);
mtf_102 = filterMassTracesForPro(mass_traces,bc,10,2);

save('mtf_51','mtf_51');
save('mtf_102','mtf_102');


[pro51,fitting51,envelopes51] = getRIAforDeconv(bc,mtf_51{1});
[pro102,fitting102,envelopes102] = getRIAforDeconv(bc,mtf_102{1});
save('filtered_results_51','pro51','fitting51','envelopes51');
save('filtered_results_102','pro102','fitting102','envelopes102');




 load('filtered_results_51.mat');
% remove all empty fits

exp_field = cellfun(@(y) cellfun(@(x) isfield(x,'expected'),y,'uni',1) , fitting51,'uni',0);

envelopes51 = cellfun(@(x,y) y(x),exp_field,envelopes51,'uni',0);
fitting51 = cellfun(@(x,y) y(x),exp_field,fitting51,'uni',0);


envelopes51 = envelopes51(~cellfun(@isempty,fitting51));
fitting51 = fitting51(~cellfun(@isempty,fitting51));

for i =1:length(pro51.reduced.peptide)

    pro51.reduced.peptide(i).ria = fitting51{i};
    pro51.reduced.peptide(i).feature = envelopes51{i};

end

peptides = pro51.reduced.peptide(cellfun(@length,fitting51) == 1);
clearvars -except peptides AAMD
for i =1:length(peptides)
    [~,atoms]=getMass_modSeq(peptides(i));
    temp = zeros(1,(atoms.C+1));
    temp(peptides(i).feature{1}.isotope + 1) = peptides(i).ria{1}.expected;
    
    peptides(i).measured = temp;
    peptides(i).weight = sum(peptides(i).feature{1}.ions(:))/sum(cellfun(@(x) sum(x{1}.ions(:)), {peptides.feature}));
end


deconv51 = automated_deconv(peptides,1,cell2mat(AAMD)');


exp_field = cellfun(@(y) cellfun(@(x) isfield(x,'expected'),y,'uni',1) , fitting102,'uni',0);

envelopes102 = cellfun(@(x,y) y(x),exp_field,envelopes102,'uni',0);
fitting102 = cellfun(@(x,y) y(x),exp_field,fitting102,'uni',0);


envelopes102 = envelopes102(~cellfun(@isempty,fitting102));
pro102.reduced.peptide = pro102.reduced.peptide(~cellfun(@isempty,fitting102));
fitting102 = fitting102(~cellfun(@isempty,fitting102));


for i =1:length(pro102.reduced.peptide)

    pro102.reduced.peptide(i).ria = fitting102{i};
    pro102.reduced.peptide(i).feature = envelopes102{i};

end

peptides102 = getPercentAgreement(pro102.reduced.peptide,deconv51.best.x);
peptides = peptides102(cellfun(@(x) sum(x>0.5) == 1, {peptides102.pa}));
peptides = reduce_features(peptides,0.5);


for i =1:length(peptides)
    [~,atoms]=getMass_modSeq(peptides(i));
    temp = zeros(1,(atoms.C+1));
    temp(peptides(i).feature{1}.isotope + 1) = peptides(i).ria{1}.expected;
    
    peptides(i).measured = temp;
    peptides(i).weight = sum(peptides(i).feature{1}.ions(:))/sum(cellfun(@(x) sum(x{1}.ions(:)), {peptides.feature}));
end

%% deconvolve at 50 percent agreement
deconv102_50 = automated_deconv(peptides,1,deconv51.best.x);
peptides102_p50 = getPercentAgreement(pro102.reduced.peptide,deconv102_50.best.x);

peptides = peptides102_p50(cellfun(@(x) sum(x>0.75) == 1, {peptides102_p50.pa}));
peptides = reduce_features(peptides,0.75);


for i =1:length(peptides)
    [~,atoms]=getMass_modSeq(peptides(i));
    temp = zeros(1,(atoms.C+1));
    temp(peptides(i).feature{1}.isotope + 1) = peptides(i).ria{1}.expected;
    
    peptides(i).measured = temp;
    peptides(i).weight = sum(peptides(i).feature{1}.ions(:))/sum(cellfun(@(x) sum(x{1}.ions(:)), {peptides.feature}));
end

%% deconolve at 75 %
deconv102_75 = automated_deconv(peptides,1,deconv102_50.best.x);


peptides102_p75 = getPercentAgreement(pro102.reduced.peptide,deconv102_75.best.x);
peptides = peptides102_p75(cellfun(@(x) sum(x>0.9) == 1, {peptides102_p75.pa}));
peptides = reduce_features(peptides,0.9);


for i =1:length(peptides)
    [~,atoms]=getMass_modSeq(peptides(i));
    temp = zeros(1,(atoms.C+1));
    temp(peptides(i).feature{1}.isotope + 1) = peptides(i).ria{1}.expected;
    
    peptides(i).measured = temp;
    peptides(i).weight = sum(peptides(i).feature{1}.ions(:))/sum(cellfun(@(x) sum(x{1}.ions(:)), {peptides.feature}));
end

%% deconolve at 90%

deconv102_90 = automated_deconv(peptides,1,deconv102_75.best.x);
peptides102_p90 = getPercentAgreement(pro102.reduced.peptide,deconv102_90.best.x);

peptides = reduce_features(peptides102_p90,0.9);

peptides = peptides(~cellfun(@isempty,{peptides.pa}));

save('final_peptides_set','peptides')