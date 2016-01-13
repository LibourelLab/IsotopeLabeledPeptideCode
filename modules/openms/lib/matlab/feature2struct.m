function [ feature ] = feature2struct( file )
%EXTRACTFEATURETXT Summary of this function goes here
%   Detailed explanation goes here

data = importdata(file);
[fids,maptofids,maptodata]=unique(data.textdata);

for i =1:length(fids)
    idx = find( i == maptodata);
    feature(i).id = fids{i};
    nr = data.data(idx,1);
    [nr_unique,map2nru,map2nrRaw]=unique(nr);
    
    for j =1:length(nr_unique)
        idx2 = idx((map2nrRaw == 1));
        mt(j).x = data.data(idx((map2nrRaw == j)),2);
        mt(j).y = data.data(idx((map2nrRaw == j)),3);
    end
    feature(i).mt = mt;
    clear mt;
end

