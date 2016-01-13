function [ out ] = isFeatureInIonRange(feature_ions,unlabaled_ions,numsd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
avg = -0.2244;
sd = 0.8582;


error = (feature_ions - unlabaled_ions)/unlabaled_ions;
max_error = avg + numsd*sd;
min_error = avg - numsd*sd;


out = error < max_error & error> min_error;


end

