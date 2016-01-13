function [ idx ] = getDendrogramDataSets(Z)
%UNTITLED9 Summary of this function goes here
%   Z is output from linkage and index is the set of the linkage function

numElements = size(Z,1)+1;

idx = cell(1,2*numElements-1);


for i = 1:(size(Z,1)+1)
    idx{i} = i; 
end

for j = numElements+1:length(idx)
   idx{j} = cell2mat(idx(cell2mat(idx(Z(j-numElements,[1:2])))));
   
    
end


end

