function [ idx ] = isMeasurementInArray(mzArray,mzMeasurement,ppm)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
idx=(cellfun(@(x) find(abs(mzArray - x)./x .*1e6 < ppm),num2cell(mzMeasurement),'uni',0));

end

