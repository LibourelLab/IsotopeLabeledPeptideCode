function [protein] = getProteinStructFromMaxQuantOutput(proteins,peptides,msms,evidence,fastafile)
% import before: 
%               proteins = dataset2struct(*import(proteins.txt); - can be
%               reduced set of interest
%               peptides =  ... cannot be reduced!
%               msms = ...     cannot be reduced!

emap = containers.Map([evidence.id],1:length(evidence));

for i =1:length(proteins)
    peptide_idx = cellfun(@str2num,regexp(proteins(i).PeptideIDs,'\;','split'),'uni',1);
    mapped_idx = cellfun(@(x) find(double([peptides(:).id]) == x),num2cell(peptide_idx),'uniform',1);
    protein(i).an = regexp(proteins(i).ProteinIDs,';','split');
    protein(i).Uniquepeptides = proteins(i).UniquePeptides;
    protein(i).Razoruniquepeptides = proteins(i).PeptideCounts_razor_unique_;
    protein(i).numPeps = proteins(i).Peptides;
    protein(i).Fastaheaders = proteins(i).FastaHeaders;
    
    %reduced peptied set to new peptide structure
    peps_i = peptides(mapped_idx);
    mapped_idx_msms = cellfun(@(x) find(double([msms(:).id]) == x),{peps_i(:).BestMS_MS},'uniform',1);
    pep_msms = msms(mapped_idx_msms);
    for j = 1:proteins(i).Peptides
         peptide(j).seq = pep_msms(j).Sequence;
         peptide(j).charge = pep_msms(j).Charge;
         peptide(j).mods = pep_msms(j).Modifications;
         peptide(j).modseq = pep_msms(j).ModifiedSequence;
         peptide(j).scannum = pep_msms(j).ScanNumber;
         peptide(j).scanidx = pep_msms(j).ScanIndex;
         peptide(j).rt  = pep_msms(j).RetentionTime;
         peptide(j).mz = pep_msms(j).m_z;
         peptide(j).mass = pep_msms(j).Mass;
         peptide(j).MassError_corr = pep_msms(j).MassError_ppm_; 
         peptide(j).MassError = pep_msms(j).SimpleMassError_ppm_;
         peptide(j).isotope = pep_msms(j).IsotopeIndex;
         peptide(j).score = pep_msms(j).Score;
         peptide(j).neutralloss = pep_msms(j).NeutralLossLevel;
         peptide(j).proteins = regexp(peps_i(j).Proteins,'\;', 'split');
         peptide(j).numProteins = length(peptide(j).proteins);
         peptide(j).isunique2group = isempty(regexp(peps_i(j).Unique_Groups_, 'no', 'once'));
         peptide(j).isunique2protein = isempty(regexp(peps_i(j).Unique_Proteins_, 'no', 'once'));
         u = evidence(emap(pep_msms(j).EvidenceID)).Uncalibrated_CalibratedM_z_ppm_;
         mz = evidence(emap(pep_msms(j).EvidenceID)).m_z;
         peptide(j).mzRef = mz*(1e-6)*u + mz;

         
    end
    protein(i).peptide = peptide;
    clear peptide;
    % map an numbers to fasta file
    protein(i).fastaMap = cellfun(@(x) find(not(cellfun(@isempty,regexp({fastafile(:).Header},x,'match','once'),'uni',1))),protein(i).an,'uni',0);
    
end






end