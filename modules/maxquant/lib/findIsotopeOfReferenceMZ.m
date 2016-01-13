function [ peptide ] = findIsotopeOfReferenceMZ(peptide)
%findIsotopeOfReferenceMZ determines the closest isotope varient to the
%reference m/z of the modified peptide sequence
%from MSMS data output.  Theoretical isotopomers are generated from the
%sequence and the closest isotopomer index is determined and appended to
%thte data structure

%input: peptide structure from maxquant output within protein data structure


mz_0 = mass2mz(getMass_modSeq(peptide),peptide.charge);
%note: atomic comp only returns number of carbons in non-modified
%sequences.  Assuming there are N carbons nonmodified sequence and N+p
%carbons in the modified sequence, we assume there wont be observable
%measurements for the last p+1 mass isotopomers.
atoms = atomiccomp(peptide.seq);
isotopes = [-3:atoms.C+1];
mzArray = mzArraySimulator('mz',mz_0,peptide.charge,1.003355,'simulate isotopes',0,isotopes);

[~,i]=min(abs(mzArray-peptide.mzRef));
peptide.isotope = isotopes(i);
end



