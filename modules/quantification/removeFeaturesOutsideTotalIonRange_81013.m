function [ clusters ] = removeFeaturesOutsideTotalIonRange_81013(clusters,protein_database,numSD)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for m=1:length(protein_database)
    
    for i =1:length(clusters{m})
    if ~isempty(clusters{m}{i})
        for j =1:length(clusters{m}{i})
            SN = featureStruct2mat(clusters{m}{i}{j},'SN',0); 
            idx = find(sum(SN,1) > 0.5*max(sum(SN,1)));
            sn=SN(:,min(idx):max(idx));
            total_ions{m}{i}{j} = sum([sn(:)]);
        end
            
            
    else
        total_ions{m}{i} = [];
    end
end

end


for m=1:length(protein_database)
    
    for i =1:length(clusters{m})
        if ~isempty(clusters{m}{i}) & ~isempty(protein_database(m).peptide(i).total_ions)
            for j =1:length(clusters{m}{i})
                q{m}{i}(j) = isFeatureInIonRange(total_ions{m}{i}{j},protein_database(m).peptide(i).total_ions,numSD);
            end
            
            clusters{m}{i} = clusters{m}{i}(q{m}{i});
            %groups{m}{i} = groups{m}{i}(q{m}{i});
            
            
        elseif isempty(protein_database(m).peptide(i).total_ions)
            q{m}{i} = NaN;
        else
 
           q{m}{i} = [];
        end
        
    end
end

end

