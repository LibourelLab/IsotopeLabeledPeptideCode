function [ mzPred,isotope_index] = simulateFullmzArray(peptide,min_iso)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

atoms = atomiccomp(peptide.seq);
isotope_index = [min_iso:atoms.C];
mzPred = mzArraySimulator('mz',peptide.mzRef,peptide.charge,1.00335,'simulate isotopes',peptide.isotope,isotope_index);
mzArray = [];
%i = peptide.charge;
%    mzArray = [mzArray,mzArraySimulator('mz',peptide.mzRef,i,1.00335,'within m/z range',min(mzPred),max(mzPred))];

end

