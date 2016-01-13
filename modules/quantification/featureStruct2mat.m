function [ M,scans] = featureStruct2mat( features,in,NonMeasurements)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mt = [features(:).mt];
SN = cellfun(@rdivide, {mt(:).intensity},{mt(:).noise},'uni',0);
MZ = {mt(:).mz};

% generate "mean" value for distrbution in feature set from the sum of S/N
% find all scans in features set
scans ={mt(:).scanNum};
scans = unique([scans{:}]);
for i =1:length(mt)
    for j = 1:length(scans)
        k = find( mt(i).scanNum == scans(j));
        if ~isempty(k)
            if ( in =='SN')
            M(i,j)=SN{i}(k);
            elseif (in =='mz')
                M(i,j)=MZ{i}(k);
            end
        else
            if NonMeasurements == false
               M(i,j) = 0;
            else
               M(i,j) = NaN;
            end
               
        end
    end
end

% for i = 1:length(mt)
%     for j = 1:length(mt)
%         L(i,j)=sum(logical(M(i,:)) & logical(M(j,:)));
%     end
% end
        
        


end

