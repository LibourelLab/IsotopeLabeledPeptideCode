
pro_red = pro;
mtf = mtf_51;


protein.reduced = pro_red; for i =1:length(protein.reduced.peptide) protein.reduced.peptide(i).idx = i; end
mass_trace_filtered.reduced = mtf;
protein.reduced.peptide(cellfun(@(x) length(x) < 2,mass_trace_filtered.reduced ,'uni',1)) = [];
mass_trace_filtered.reduced(cellfun(@(x) length(x) < 2,mass_trace_filtered.reduced ,'uni',1)) = [];
    
[S]= cellfun(@(x,y) getSimilarityMatrix_80713([x.mt],protein.reduced.peptide(y),4,1),mass_trace_filtered.reduced,num2cell([1:length(protein.reduced.peptide)]), 'uni',0);
Z = cellfun(@(x) linkage(x,'ward','euclidean'),S,'uni',0);
groups = cellfun(@(x,y) clusterMassTracesFromDendrogram(x,[y.mt],2.5),Z,mass_trace_filtered.reduced,'uni',0);
mt = cellfun(@(x) [x.mt],mass_trace_filtered.reduced,'uni',0);
%sort group

x = cellfun(@isempty,groups, 'uni',1);
mass_trace_filtered.reduced(x) = []; protein.reduced.peptide(x) = []; mt(x) = []; groups(x) = [];


groups = cellfun(@(x) cellfun(@sortGroup,x,'uni',0),groups,'uni',0);

clusters = cellfun(@(y,z) cellfun(@(x) y(x.sets),z,'uni',0),mt,groups,'uni',0);

%break clusters into subfeatures with 
clusters = cellfun(@(x) breakCluster(x,2,1),clusters,'uni',0);
x = cellfun(@length,clusters,'uni',1) == 0;
mass_trace_filtered.reduced(x) = []; protein.reduced.peptide(x) = []; mt(x) = []; groups(x) = []; clusters(x) = [];

clust_pro{1} = clusters;


% [ clust_pro ] = removeFeaturesOutsideTotalIonRange_81013(clust_pro,protein.reduced,2);
% clusters = clust_pro{1}; 

[clusters,n]=cellfun(@removeClustersWithNegativeIsotopes_81013,clusters,'uni',0);
[clusters,n]=cellfun(@removeFeaturesOfDiffChargeState_81013,clusters,'uni',0);

x = cellfun(@isempty,clusters,'uni',1);
clusters(x) = [];
protein.reduced.peptide(x) = [];



envelopes = cellfun(@(z,y) cellfun(@(x) getEnvelope(x,protein.reduced.peptide(y),60000,60000,6,0.5,40),z,'uni',0),clusters,num2cell([1:length(clusters)]),'uni',0);



%fitting_red = fitting;
envelopes_red = envelopes;

% for i =1:length(envelopes_red)
for i =length(clusters)
     k = ~logical(cellfun(@(x) sum([x.ions(:)]) || length(x.isotope)<2, envelopes_red{i},'uni',1));
    envelopes_red{i}(k) = [];
end
fitting = cell(1,length(envelopes_red));
for i =1:length(envelopes_red)
    for j=1:length(envelopes_red{i})
        
       fitting{i}{j} = fitRIA2ExperimentalData(envelopes_red{i}{j}.ions);
    end
end


for i =1:length(envelopes_red)
   k =  cellfun(@(x) length(x.isotope),envelopes_red{i},'uni',1) == 1;
   envelopes_red{i}(k) = [];
   fitting{i}(k) = [];
%     residuum_l1{i} = getFeatureResiduum_81113(fitting{i},envelopes_red{i},protein.reduced.peptide(i));
     residuum_l1{i} = getFeatureResiduum_81113(fitting{i},envelopes_red{i},pro_1223.peptide(i));

end





for i =1:length(envelopes)
    
    for j=1:length(envelopes{i})
        if ~isempty(envelopes{i}{j}.ions)

         
         if ~isempty(fitting{i}{j}) & ~negative_isos{i}(j)
             if length(envelopes{i}{j}.isotope) > 1
             isoMAP = envelopes{i}{j}.isotope;
             
             meas_mid = zeros(1,length(sev.peptide(i).simulatedMID));
             meas_mid(isoMAP+1) = fitting{i}{j}.expected;
             
            residuum{i}{j} = sum((sev.peptide(i).simulatedMID - meas_mid).^2);
             else
                   residuum{i}{j} = [];
             end
         else
             residuum{i}{j} = [];
         end
        else
                         residuum{i}{j} = [];
        end
    end
    
end

string = strcat('mass_weight_r_',num2str(n));
save(string,'S','groups','clusters','envelopes','fitting','residuum');

end

    

