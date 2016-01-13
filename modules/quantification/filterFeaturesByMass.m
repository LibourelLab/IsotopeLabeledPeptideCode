function [ features, idx] = filterFeaturesByMass(features,mzArray,ppm)
%UNTITLED11 Summary of this function goes here

% 1st, remove all features not within the specied range of mzArray
    mzMIN = min(mzArray) - min(mzArray)*ppm/1e6;
    mzMAX = max(mzArray) + max(mzArray)*ppm/1e6;


    mzRanges = cellfun(@(x) [min(x), max(x)], {features(:).mz},'uni',0);
    inMZrange = cellfun(@(x) x(2) > mzMIN & x(1) < mzMAX, mzRanges,'uni',1);
    features= features(inMZrange);
    idx = find(inMZrange);

end

