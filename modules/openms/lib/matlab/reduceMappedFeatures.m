function [ features_mapped ] = reduceMappedFeatures(features_mapped)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

for i =1:length(features_mapped)
   for j =1:length(features_mapped(i).mt)
       weights = features_mapped(i).mt(j).intensity ./ sum(features_mapped(i).mt(j).intensity);
       features_mapped(i).mz(j) = weights*features_mapped(i).mt(j).mz';
       features_mapped(i).rt(j) = weights*features_mapped(i).mt(j).rt';
   end 



end

end