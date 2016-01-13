function [ residual ] = fitExpectedValue2MeasuredIons(x,A)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%make sure x is horizontal vector
[a1,~]=size(x);
l = length(x);

if (a1 == l)
    x = x';
else
end
    


%important that A is arranged such that rows correspond to mass and column
%correspond to retention time

[numMasses,numScans] = size(A);
Mask = double(logical(A));
% parition x into expected value vector, E and number of ions per scan, N
E = x(1:numMasses);
N = x(numMasses+1: numMasses+numScans);

residual = 0;
for i = 1:numScans
    
    %residual = residual + sum((N(i)*E' - A(:,i)).^2.*Mask(:,i));
    %residual = residual + sum((N(i)*E' - A(:,i)).^2);
    residual = residual + sum(abs(N(i)*E' - A(:,i)).*Mask(:,i));
end


end

