function [Z,z_sets,r] = reduceDendrogram(Z,mass_traces)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

z_sets = getDendrogramDataSets_8313(Z,length(mass_traces));
z_sets([1:length(mass_traces)]) = [];

for i =1:length(z_sets)
    isotopes = [mass_traces(z_sets{i}).isoIdx];
    if length(isotopes)  == length(unique(isotopes))
        req(i) = true;
    else
        req(i) = false;
    end


end

Z=Z(req,:);
z_sets = z_sets(req);
r = find(req);


end

