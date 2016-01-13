function [S] = getSimilarityMatrix(mass_traces,peptide,ppm_cutoff,r)

%full mass traces structure
%s/n weighted m/z value for entire mass traces
if length(mass_traces) > 1 
mz = [mass_traces.mz];

[SN,scans] = featureStruct2mat(mass_traces,'SN',0);
[MZ] = featureStruct2mat(mass_traces,'mz',0);

C = cov(SN');
SD = sqrt(diag(C));
R = C./(SD*SD');
R = (R-min(R(:))) ./ (max(R(:)-min(R(:))));




[~,~,ppm] = mzArrayIsotopicDistance(mz,peptide.mzRef,peptide.charge,1.003355);
u = ppm_cutoff*sqrt(1/(2*log(2)));
W = gaussmf(squareform(pdist(ppm','cityblock')),u,0);

%S=getSimilarityMatrix(R,W);
S = rescaleMatrix(W,r).*R;

elseif length(mass_traces) == 0
    S= [];
else
    S = 1;
end
% 
% data.S = S;
% data.corr = R;
% data.cov = C;
% data.weight = W;
% data.ppm = ppm;
% data.SN = SN;
% data.scans = scans;
% data.mz = mz;
% data.MZ = MZ;


end

