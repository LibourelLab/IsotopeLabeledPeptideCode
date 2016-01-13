function [ group ] = sortGroup(group)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,k]=sort([group.sets]);

group.sets = group.sets(k);
end

