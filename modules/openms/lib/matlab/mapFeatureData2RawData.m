function [ feature_out ] = mapFeatureData2RawData(rawData,feature,isNoiseInRawData)
feature_out = struct('mt',[]);
mt = struct('scanNum',[],'rt',[],'idx',[],'mz',[],'intensity',[],'noise',[]);

for i =1:length(feature)
   for j = 1:length(feature(i).mt)
       for k = 1:length(feature(i).mt(j).x)
             rt = feature(i).mt(j).x(k);
             mz =  feature(i).mt(j).y(k);
                    [~,q]= min(abs(rawData.RT.*60 - rt));
                    scan = rawData.scans(q);
                    [~,p] = min(abs(scan.mz - mz));
                    mt(j).scanNum(k) = q;
                    mt(j).rt(k) = rawData.RT(q);
                    mt(j).idx(k) = p;
                    mt(j).mz(k) = scan.mz(p);
                    mt(j).intensity(k) = scan.intensity(p);
                    if(isNoiseInRawData)
                        mt(j).noise(k) = scan.noise(p);
                    end
                    
       end
       
       [~,m] =unique([mt(j).scanNum]);
       mt(j).scanNum = mt(j).scanNum(m);
       mt(j).rt =  mt(j).rt(m);
       mt(j).idx =  mt(j).idx(m);
       mt(j).mz = mt(j).mz(m);
       mt(j).intensity = mt(j).intensity(m);
       if(isNoiseInRawData)
          mt(j).noise = mt(j).noise(m);
       end
   end
   feature_out(i).mt = mt;
   clear mt;
end


end

