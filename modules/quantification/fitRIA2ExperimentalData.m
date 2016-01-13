function [ out ] = fitRIA2ExperimentalData(IonCountMatrix)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
out = [];
if ~isempty(IonCountMatrix)
    

[numMasses,numScans] = size(IonCountMatrix);
if numScans > 1 & numMasses > 1
expected = sum(IonCountMatrix')./sum(sum(IonCountMatrix));
numIons_scan =  sum(IonCountMatrix);

x0 = [expected,numIons_scan];

Aeq = [ones(1,numMasses),zeros(1,numScans)];
beq = 1;
lb = zeros(1,numMasses+numScans);
ub = [ones(1,numMasses),inf*ones(1,numScans)];
A = -[zeros(1,numMasses),ones(1,numScans)];
b = -sum(sum(IonCountMatrix));

% Sets options for optimizer
options = optimset('Display','Iter-detailed','TolFun',1e-11,'TolX',1e-11, ...
                    'MaxFunEval',1e30,'MaxIter',1e4,'Algorithm','interior-point','DiffMinChange',1e-4);


[x_opt,fval] = fmincon(@fitExpectedValue2MeasuredIons,x0,A,b,Aeq,beq,lb,ub,[],options,... 
        IonCountMatrix);
    
    
out.expected = x_opt(1:numMasses).*double(logical(sum(IonCountMatrix')));
out.numIons_scan = x_opt(numMasses+1:numMasses+numScans);
out.fval = fval;
    
elseif numMasses == 1
    out = [];
    
else
    out.expected = IonCountMatrix./sum(IonCountMatrix);
    out.numIons_scan = sum(IonCountMatrix);
    out.fval = 0;
    
end
end


end

