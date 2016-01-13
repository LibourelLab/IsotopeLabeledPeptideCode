function [ cluster ] = sortCluster(cluster)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,k]=sort([cluster.isoIdx]);
cluster = cluster(k);

end

