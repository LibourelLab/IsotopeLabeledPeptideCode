function [peptide] = cluster_mass_traces(pro,mtf)

% allocated protein structure to new protein.reduced structure (useful for
% troubleshooting)
protein.reduced = pro; 

% for i each peptide add an index
for i =1:length(protein.reduced.peptide) 
    protein.reduced.peptide(i).idx = i; 
end

% make a similar mass trace filtered structure
mass_trace_filtered.reduced = mtf;


% remove all peptides with less than 2 mass traces (no isotope distribution
% information is contained within a single peak)
protein.reduced.peptide(cellfun(@(x) length(x) < 2,mass_trace_filtered.reduced ,'uni',1)) = [];

% repeat previous step on mass trace structure
mass_trace_filtered.reduced(cellfun(@(x) length(x) < 2,mass_trace_filtered.reduced ,'uni',1)) = [];
    
% get Similarity matrix for each set of mass traces (to be use for
% heiarchical clustering)
S = cellfun(@(x,y) getSimilarityMatrix([x.mt],protein.reduced.peptide(y),4,1), ...
    mass_trace_filtered.reduced,num2cell(1:length(protein.reduced.peptide)), 'uni',0);

% get collection of dendrograms
Z = cellfun(@(x) linkage(x,'ward','euclidean'),S,'uni',0);

% call clustering algorithm from dendrogram, using the emperically derived
% tree cutoff of 2.5 (see MS thesis, Automated Identification of 13C labeled
% labeled peptides for description of paramter choice)

groups = cellfun(@(x,y) clusterMassTracesFromDendrogram(x,[y.mt],2.5),Z,mass_trace_filtered.reduced,'uni',0);

% create simpler mass trace group for filtering in subsequent steps
mt = cellfun(@(x) [x.mt],mass_trace_filtered.reduced,'uni',0);


% find a remove all peptides with no cluster
x = cellfun(@isempty,groups, 'uni',1);
mass_trace_filtered.reduced(x) = []; protein.reduced.peptide(x) = []; mt(x) = []; groups(x) = [];

% sort groups of mass traces
groups = cellfun(@(x) cellfun(@sortGroup,x,'uni',0),groups,'uni',0);

% create clusters from groups of mass traces
clusters = cellfun(@(y,z) cellfun(@(x) y(x.sets),z,'uni',0),mt,groups,'uni',0);

%break clusters into subfeatures based on gapsize filter (gap-size
%paramteer = 2 (see MS thesis, Automated Identification of 13C labeled
% labeled peptides for description of paramter choice)
clusters = cellfun(@(x) breakCluster(x,2,1),clusters,'uni',0);

% find peptides with no clusters after gapsize filter
x = cellfun(@length,clusters,'uni',1) == 0;
mass_trace_filtered.reduced(x) = []; protein.reduced.peptide(x) = []; mt(x) = []; groups(x) = []; clusters(x) = [];


% if the total ion number was counted, remove features with wildly
% different ion counts between reference and labeled samples
if isfield(protein.reduced,'total_ions')
    clust_pro{1} = clusters;
    [ clust_pro ] = removeFeaturesOutsideTotalIonRange_81013(clust_pro,protein.reduced,2);
    clusters = clust_pro{1};
end


% remove features with negative isotopes
[clusters,n]=cellfun(@removeClustersWithNegativeIsotopes_81013,clusters,'uni',0);

% remove features with different charges states
[clusters,n]=cellfun(@removeFeaturesOfDiffChargeState_81013,clusters,'uni',0);


% remove empty clusters,peptides
x = cellfun(@isempty,clusters,'uni',1);
clusters(x) = [];
protein.reduced.peptide(x) = [];


% get envelopes (ion count matrices)
envelopes = cellfun(@(z,y) cellfun(@(x) getEnvelope(x,protein.reduced.peptide(y),60000,60000,6,0.5,40),...
    z,'uni',0),clusters,num2cell([1:length(clusters)]),'uni',0);




% for each peptide, eliminate envelopes with < 2 mass traces
for i =length(clusters)
     k = ~logical(cellfun(@(x) sum([x.ions(:)]) || length(x.isotope)<2, envelopes{i},'uni',1));
    envelopes{i}(k) = [];
end

% create an array for all fitted relative isotope abundances
fitting = cell(1,length(envelopes));

% for each peptide
for i =1:length(envelopes)
    % for each candidate feature
    for j=1:length(envelopes{i})
       % fit releative isotope to experimental data
       fitting{i}{j} = fitRIA2ExperimentalData(envelopes{i}{j}.ions);
    end
end


%% create an output structure for each protein/envelope/RIA fit
peptide = protein.reduced.peptide;


for i = 1:length(peptide)
    
    % find envelopes which no expectation could be formulated from fit ria
    k = find(cellfun(@isempty,envelopes{i}) | cellfun(@isempty,fitting{i}));
    
    % get rid of empty containers
    envelopes{i}(k) = [];    
    fitting{i}(k) = [];
end

for i = 1:length(peptide)
    peptide(i).feature = struct();
    for j = 1:length(envelopes{i})
       
        peptide(i).feature(j).ria = fitting{i}{j}.expected;
        peptide(i).feature(j).numIons_scan = fitting{i}{j}.numIons_scan;
        peptide(i).feature(j).fval = fitting{i}{j}.fval;
        peptide(i).feature(j).ions = envelopes{i}{j}.ions;
        peptide(i).feature(j).scans = envelopes{i}{j}.scans;
        peptide(i).feature(j).mz = envelopes{i}{j}.mz;
        peptide(i).feature(j).isotope = envelopes{i}{j}.isotope;
        peptide(i).feature(j).mz_mat = envelopes{i}{j}.mz_mat;
        peptide(i).feature(j).fid = envelopes{i}{j}.fid;
        peptide(i).feature(j).ppm = envelopes{i}{j}.ppm;

    end


end
