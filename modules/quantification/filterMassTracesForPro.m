function [ mass_traces_filtered ] = filterMassTracesForPro(mass_trace,protein,MAC,RTW)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
mass_traces_filtered = cell(1,length(protein));
for i =1:length(protein)
    mass_traces_filtered{i} = cell(1,length(protein(i).peptide));
    for j =1:length(protein(i).peptide)
        mass_traces_filtered{i}{j} = filterMassTracesForPeptide(mass_trace,protein(i).peptide(j),MAC,RTW);
    end
end

end

