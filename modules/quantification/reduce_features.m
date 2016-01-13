function [peptides] = reduce_features(peptides,cutoff)
%% reduce features based on percent agreement cutoff
% JEG 11-12-13


% for each peptide
for i =1:length(peptides)
    % construct a vector indicating which features are going to be kept
    keep = true(length(peptides(i).feature),1);
    % for each feature
    for j = 1:length(peptides(i).feature)
        % if percent agrement for candidate features is below cutoff
        if peptides(i).feature(j).percent_agreement < cutoff
            % record for deletion later
            keep(j) = false;
        end
    end
    % only keep candidate features above percent agreement threshold 
    peptides(i).feature = peptides(i).feature(keep);
end
end