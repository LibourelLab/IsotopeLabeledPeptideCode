function [ groupedInt ] = groupIntegers(array, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(varargin)
    gapsize = 1;
else
    gapsize = varargin{1};
end
    
groupedInt = {[]};
if (~isempty(array))
groupedInt = cell(1);
gIdx =1;
k = array(1);
groupedInt{gIdx} = k;

for i =2:length(array)
   if (array(i) <= k + gapsize) 
       groupedInt{gIdx}(end+1) = array(i);
   else
       gIdx = gIdx + 1;
       groupedInt{gIdx} = array(i);
   end
   k = array(i);
   
end
else

end

