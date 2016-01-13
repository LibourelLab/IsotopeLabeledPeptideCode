function [ out ] = breakCluster(clusters,gapsize,removeSize1)
% cluster cell structure: cluster{i} for ith protein, cluster{i}{j} for jth
% peptide of protein i, cluster {i}{j}{k} == kth cluster
if(~isempty(clusters))

clusters = cellfun(@sortCluster,clusters,'uni',0);
out = cell(1,length(clusters));

        for k = 1:length(clusters)
            
           if iscell(clusters)
                
           g = groupIntegers([clusters{k}.isoIdx],gapsize);
           clusters{k} = mat2cell(clusters{k},1,cellfun(@length,g,'uni',1));
           out = [out, clusters{k}];
            else
                
           g = groupIntegers([clusters.isoIdx],gapsize);
           clusters = mat2cell(clusters,1,cellfun(@length,g,'uni',1));
           out = [out, clusters];
                
        end
        out(1) = [];
        
        if removeSize1
            out(cellfun(@length,out,'uni',1) == 1) = [];
        end
        end
else
    out = [];
end
        
end

