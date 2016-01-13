function [ filtered_mass_trace ] = filterMassTraceForProtein(mass_trace,protein,accuracy,rt_range)
filtered_mass_trace = cell(1,length(protein));

for i = 1:length(protein)
    filtered_mass_trace{i} = cell(1,length(protein(i).peptide));
    
    for j = 1:length(protein(i).peptide)
        filtered_mass_trace{i}{j} = filterMassTracesForPeptide_71313(mass_trace,protein(i).peptide(j),accuracy,rt_range);    
        
    end
end



end

