function [isotope ] = joinBrokenMassTraces(features,ppmCutoff,scanBreak)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%group mass traces by isotope index

[isos,u]=sort(([features(:).isoIdx]));
isos_u = unique(isos);

idx=mat2cell(u',cellfun(@(y) sum(isos == y), num2cell(isos_u),'uni',1)',1);
isotope = struct('mt',[]);
for k =1:length(idx)
    mt = features(idx{k});
    [~,s]=sort(cellfun(@min, cellfun(@(x) x.scanNum,{mt(:).mt},'uni',0),'uni',1));
    mt = mt(s);
    i = 1;
    while (i < length (mt))
        j = i+1;
        while (j <= length(mt))
            c0 = abs(mt(i).mz - mt(j).mz)./mt(i).mz.*1e6 < ppmCutoff;
            c1 = mt(i).mt.scanNum(end) < mt(j).mt.scanNum(1);
            c2 = mt(j).mt.scanNum(1)- mt(i).mt.scanNum(end) < scanBreak;
            if(c1&c2&c0)
                mt(i).mt = joinMassTraces(mt(i).mt,mt(j).mt);
                mt(j) = [];
            else
                j = j+1;
            end
        end
    i = i+1;
    end
    
    for i =1:length(mt)
       
        mt(i).mz = sum(mt(i).mt.mz.*mt(i).mt.intensity./mt(i).mt.noise)./sum(mt(i).mt.intensity./mt(i).mt.noise);
        mt(i).rt = sum(mt(i).mt.rt.*mt(i).mt.intensity./mt(i).mt.noise)./sum(mt(i).mt.intensity./mt(i).mt.noise);
        
    end
    
    isotope(k).mt = mt;
    
    
    
end



end