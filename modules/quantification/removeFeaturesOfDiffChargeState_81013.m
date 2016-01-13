function [ clusters,k ] = removeFeaturesOfDiffChargeState_81013(clusters)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if(~isempty(clusters))
    
for i =1:length(clusters)
    
   u2 = unique(mod([clusters{i}.isoIdx],2));
   u3 = unique(mod([clusters{i}.isoIdx],3));

   if length(u2)<2 || length(u3) <3
       k(i) = true;
   else
       k(i) = false;
   end
   
    
end
clusters(k) = [];
else
    k = [];
    clusters = [];

end

