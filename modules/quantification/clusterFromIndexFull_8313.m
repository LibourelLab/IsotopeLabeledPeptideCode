function [ out ] = clusterFromIndexFull_8313(Z_i,Z_f,mapping,idx,numNodes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Z is dendrogram coordinates
% corr is the group metric used to define inner group variance, currently
% just mean correlation
% idx is the group index we wish to group from
% cutoff is the minumum correlation a group should have

%numNodes = size(Z_i,1)+1; 

%group all node cooridates
Z_c=num2cell(Z_i(:,1:2),2);

%get index mapping at each node (i.e. all data elements (index sets) contained at each
%node, 
%Z_f=getDendrogramDataSets(Z_i);
%Z_f(1:numNodes) = [];


k = 0;
x = true;
while(x)
k = k+1;
%find node that contains index and record contents in idx_f{k}
temp=find(cellfun(@(x) any(x == idx),Z_c,'uni',1));
if ~isempty(temp)
idx_i(k)=find(cellfun(@(x) any(x == idx),Z_c,'uni',1));
height{k} = Z_i(idx_i(k),3);
    if idx_i(k) == length(Z_c)
        x = false;    
    else
        idx = mapping(idx_i(k)) + numNodes;
    end
else
    x = false;
end

out= struct('sets',Z_f(idx_i),'Zidx',num2cell(idx_i),'height',height);

end

    
