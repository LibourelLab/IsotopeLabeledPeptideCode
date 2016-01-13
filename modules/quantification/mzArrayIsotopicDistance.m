function [ dist,idx,signed_distance] = mzArrayIsotopicDistance_51514(mzList,mzRef,charge,cmd)
%UNTITLED2 Assuming all values in mzList are potential parts of an isotopic
%cluster inluding a reference m/z value and charge(mzRef & charge), we can 
%determine an approximate precision (in ppm) of m/z charge values within
%that list

% ----- Inputs: mzList: mass/charge array
%               mzRef: reference mass/charge
%               charge: charge of isotopic cluster
%               cmd: mass difference between carbon 12 and 13

% This function takes the modulus of the mzList and msRef with respect to
% cmd/charge.  Theses moduli are then converted to angles (in degrees).
% The distance between the angles are calculated and converted to an
% appoxiated ppm.

% NOTE: this is an approximate ppm precision because the acutual formula
% would require the true m/z value for the 
% 

% sim
dist = [];
if (isempty(cmd))
    cmd = 1.003355;
end



predMZ=mzArraySimulator('mz',mzRef,charge,cmd,'within m/z range',min(mzList),max(mzList));



%match predicted mass/charge value for each observed m/z List
for i =1:size(mzList,1)
    for j =1:size(mzList,2)

    [~,idx(i,j)]=min(abs(mzList(i,j) - predMZ));
    pred_mz_array(i,j) = predMZ(idx(i,j));
    end
end

for i =1:size(mzList,1)
    for j =1:size(mzList,2)

   if (mzList(i,j) - predMZ(idx(i)))<0
       sign(i,j) = -1;
   else
       sign(i,j) = 1;
   end
    end
end
%idx = idx - 1;



%cmd = 1.00335;
angles=mod(mzList,cmd/charge)*(charge/cmd)*360;
refAngle = mod(mzRef,cmd/charge)*(charge/cmd)*360;


for i =1:size(mzList,1)
    for j =1:size(mzList,2)
        
    [dist_i] = distanceBetweenAngles(refAngle,angles(i));
    dist(i,j) = (dist_i./(360*charge/cmd)/pred_mz_array(i))*1e6;
    %s(i) = (s_i./(360*charge/cmd)/pred_mz_array(i))*1e6
    end
end


signed_distance = sign.*dist;


end

