function [ isotopes ] = filterMassTracesForPeptide_71313(features,peptide,ppm,rt_delta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%rt_delta = 1;
%ppm = 20;
ppm_patch = 2;

%parameters for second reduction
numDiscontinuousScans = 0;
cutoff = 1.8;
maxNumDiscontinuousIsotopes = 2;



for i = 1:length(features)
    for j = 1:length(features(i).mt)
    features(i).mt(j).fid = i;
    end
end
min_iso = -3;
[mzPred,isotope_index]=simulateFullmzArray(peptide,min_iso);
mzRef = peptide.mzRef;
%first filter all features by m/z predicted range and then by retention
%time range
features=filterFeaturesByMass(features,mzPred,ppm);
features=filterFeaturesByRetentionTime(features,peptide.rt-rt_delta,peptide.rt+rt_delta);
    

%look for masses in each feature to fall within specified cutoff (i.e.
%filter all masses that might be within the predicted m/z range, but are
%not in within the mass accuracy of predicted isotopes

matched_isotopes = cellfun(@cell2mat,cellfun(@(x) isMeasurementInArray(mzPred,x,ppm),{features(:).mz},'uni',0)','uni',0);

idx = find(not(cellfun(@isempty,matched_isotopes,'uni',1)));
features=features(idx);
for i =1:length(features)
   features(i).isoIdx = isotope_index(matched_isotopes{idx(i)});
end

% join broken mass traces and group into possible isotopes
if ~isempty(features)
isotopes = joinBrokenMassTraces(features,ppm_patch,numDiscontinuousScans);
else 
    isotopes = [];
end




end

