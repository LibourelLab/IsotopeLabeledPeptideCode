function [ out ] = getEnvelope(mass_traces,peptide,R,R0,K,r,ppm_refilt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if isempty(r)
    r = 0.5;
end


out = struct('ions',[]);
charge = peptide.charge;
[SN,scans] = featureStruct2mat(mass_traces,'SN',0);
[MZ,scans] = featureStruct2mat(mass_traces,'mz',0);



idx=find(sum(SN) > r*max(sum(SN)));
SN = SN(:,min(idx):max(idx));
MZ = MZ(:,min(idx):max(idx));

IONS = SN.* K * sqrt(R0 / R) / charge;
out.ions = IONS;
out.scans = scans(min(idx):max(idx));
out.mz = [mass_traces(:).mz];
out.isotope = [mass_traces(:).isoIdx];
out.mz_mat = MZ;
out.fid = cellfun(@(x) x.fid, {mass_traces.mt});

filterd = out;
mzRef = peptide.mzRef;

[numIons,numScans] = size(MZ);
for i =1:numIons
   [~,~,PPM(i,:)] = mzArrayIsotopicDistance(MZ(i,:),peptide.mzRef,peptide.charge,1.003355);
end
%PPM = PPM(:,min(idx):max(idx));
out.ppm = PPM.*double(logical(SN));

%rescale to median ppm in cluster
p = [out.ppm(:)];
%out.ppm = out.ppm - median(p);

if ~isempty(ppm_refilt)
mask = double(abs(out.ppm) < ppm_refilt);
out.ions=out.ions.*mask;
end

i = logical(sum(out.ions'));
j = logical(sum(out.ions));
out.ions = out.ions(i,j);
out.mz = out.mz(i);
out.ppm = out.ppm(i,j);
out.isotope = out.isotope(i);
out.scans = out.scans(j);
out.mz_mat = out.mz_mat(i,j);
out.fid = out.fid(i);

end

