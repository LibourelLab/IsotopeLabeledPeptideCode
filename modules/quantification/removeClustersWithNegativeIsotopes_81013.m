function [ clusters,n] = removeClustersWithNegativeIsotopes_81013( clusters)
%UNTITLED4 Summary of this function goes c
%   Detailed explanation goes here
n = [];
if (~isempty(clusters))
n = cellfun(@(x) any([x.isoIdx]<0) ,clusters,'uni',1);
clusters(n) = [];
end


end

