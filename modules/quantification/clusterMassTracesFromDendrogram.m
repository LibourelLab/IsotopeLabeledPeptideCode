function [groups] = clusterMassTracesFromDendrogram(Z,mtf,heightcutoff)

z_spent = [];
j = 1;
sets = [];
groups = [];


[Z_i,Z_sets,mapping] = reduceDendrogram(Z,mtf);
Z_f = getDendrogramDataSets(Z);
Z_f(1:length(mtf)) = [];
m = 1:length(mtf);

while(~isempty(Z_sets))


%identify group with least variance
[~,h]= min(Z_i(:,3));
idx = Z_i(h,1);

out = clusterFromIndexFull_8313(Z_i,Z_sets,mapping,idx,length(mtf));
condition2 = find((diff([out.height]) > heightcutoff));

if ~isempty(condition2)
    k = min(condition2) + 1;
else
    k = length(out);
end

groups{j} = out(k);

%find all dendrogram index nodes that are contained within the parent node
%in groups{j}
% z_i_all = cell2mat(cellfun(@(x) Z(x,[1,2]),{out(1:k).Zidx} ,'uni',0));
% z_spent = z_i_all(z_i_all > length(mtf)) - length(mtf);
% z_spent = [z_spent,[out(1:k).Zidx]];
% z_spent = unique(z_spent);

z_spent = find(~cellfun(@(x) isempty(intersect(x,[out(1:k).sets])),Z_sets,'uni',1));

%sets = [sets,groups{j}.sets];
%z_spent =[z_spent, out(1:k).Zidx];
Z_i(z_spent,:) = [];
Z_sets(z_spent) = [];
mapping(z_spent) = [];
j = j+1;
end


if (~isempty(groups))
m = setdiff(1:length(mtf),cell2mat(cellfun(@(x) [x.sets],groups,'uni',0)));
for i = 1:length(m)
    groups{j} = struct('sets',m(i),'Zidx',[],'height',[]);
    j = j+1;
end
end



end

