function [ features, idx] = filterFeaturesByRetentionTime( features, rtMIN,rtMAX )
%UNTITLED12 Summary of this function goes f
%   Detailed explanation goes here
    
rtRanges = cellfun(@(x) [min(x), max(x)], {features(:).rt},'uni',0);
inRTrange= cellfun(@(x) x(2) > rtMIN & x(1) < rtMAX, rtRanges,'uni',1);
features= features(inRTrange);
idx = find(inRTrange);
end

