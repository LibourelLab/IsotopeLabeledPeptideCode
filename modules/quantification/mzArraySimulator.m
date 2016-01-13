function [ massArray ] = mzArraySimulator(inputType,refMass,charge,delta,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% varargin:  1st input needs to be a flag: 
%                                          - 'simulate isotopes'
%
%                                            *Mode takes two inputs:
%                                             refIndex (0 for 'm+0) and number array for
%                                             isotopes requrested ([1,2]
%                                             for 'm+1' and 'm+2'
%
%                                          - 'within m/z range'
%                                            *Mode takes two inputs:
%                                             minimum & maximum mass or m/z
% Notes: all mass inputs are assumed to NOT include proper number of
% protons for a given charge state


% if no carbon mass difference is specfied, use the cited mass difference

%constants: proton mass and carbon isotope mass difference
p_mass = 1.00727638;
massArray = [];
if(isempty(delta))
    delta = 1.003355;
end

allowableTypes = {'mz','mass/charge', 'm/z','mass','amu'};
typeMatch=regexp(inputType,allowableTypes,'ONCE');

if (isempty(typeMatch))
   error('input #2 not idetified as a mass or mass/charge value,  Please enter m/z or mass.');
end

modeMatch=regexp(varargin{1},{'simulate isotopes','within m/z range'},'ONCE');

if( isempty(modeMatch))
    error('input #5 needs to be "simulate isotopes" or "within m/z range"');
end

% find mode index (1 for simulate isotopes, 2 for within m/z range)
mode = find(not(cellfun(@isempty,modeMatch,'uniform',1)));
type = find(not(cellfun(@isempty,typeMatch,'uniform',1)));

if (type >3 )
    %if given a mass input, convert to m/z
    refMass = (refMass + p_mass*charge)/charge;        
    
end


if (mode == 1)
    [refIdx,n]=deal(varargin{2},varargin{3});
        for i=1:length(n)
            massArray(i) = refMass + (n(i)-refIdx)*delta/charge;
        end 
    

elseif (mode == 2)
    [minMZ,maxMZ]=deal(varargin{2},varargin{3});
    if (type>3)
        minMZ = (minMZ + p_mass*charge)/charge;
        maxMZ = (maxMZ + p_mass*charge)/charge;
    end
           refModulus = mod(refMass,delta/charge);
           minIdx=  floor((minMZ - refMass)/(delta/charge));
           maxIdx = ceil((maxMZ - refMass)/(delta/charge));
           range = [minIdx:maxIdx];
           
           for i = 1:length(range)
               massArray(i) = refMass + (range(i))*delta/charge;
           end
end
  
if (type > 3)
    
    massArray = massArray.*charge;
else
end
           
           
        
end








% function [mzList] = getMZlist(m0,z,n,carbonMassdiff)
% %UNTITLED4 Summary of this function goes here
% %   Detailed explanation goes here
% % m0 == base MASS
% % z == charge state
% % n == index or index set for simulated mz array
% % 
% p_mass = 1.00727638;
% %carbonMassdiff = 1.00335;
% 
% for i=1:length(n)
%       mzList(i) = m0 + n(i)*carbonMassdiff;
%         
% end
% mzList = (mzList + z*p_mass)./z;
% 
% 
% 
% 
% %% old way of calculating m/z array
% % % neutron mass
% % n_mass =  1.0086649156;
% % %proton mass
% % p_mass = 1.00727638;
% % 
% % 
% % for i=1:length(n)
% %         mzList(i) = m0 + n(i)*n_mass;
% %         
% %     end
% % mzList = (mzList + z*p_mass)./z;
% end
