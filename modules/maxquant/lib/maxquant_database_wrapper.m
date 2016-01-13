function [ protein_database ] = maxquant_database_wrapper(protein_database,mzxmldata)
%maxquant_database_wrapper: curates raw output from maxquant
%   Detailed explanation goes here

%1) check for missed modifications and append RT field
for i =1:length(protein_database)
    for j =1:protein_database(i).numPeps
    protein_database(i).peptide(j) = checkModifiedSequences(protein_database(i).peptide(j));
    protein_database(i).peptide(j).RT = protein_database(i).peptide(j).rt;
    end
end        

% re-map precursor ion and isotope index to peptide
protein_database = getPrecursorMass(protein_database,mzxmldata);

%remove all peptides with negative isotope indexs (could be contamination)
for i =1:length(protein_database)
    k = [];
    for j =1:protein_database(i).numPeps
        if protein_database(i).peptide(j).isotope < 0;
            k = [k,j];
        end
      
    end
    
    protein_database(i).peptide(k) = [];
    protein_database(i).numPeps = protein_database(i).numPeps - length(k);
end
        




end

