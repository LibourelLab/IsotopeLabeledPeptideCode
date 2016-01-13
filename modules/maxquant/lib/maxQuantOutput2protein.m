function [ protein ] = maxQuantOutput2protein(directory,fastafile)
%maxQuant_parser: Function that imports data from maxquant text directory and fasta file to generate protein data structure 
%   Inputs: directory = working directory with all maxQuant ouutput
%           fasta - fasta file used in proteome search
fasta = fastaread(fastafile);
msmsfile = strcat(directory,'msms.txt');
peptidesfile = strcat(directory,'peptides.txt');
proteinGroupsfile = strcat(directory,'proteingroups.txt');
evidencefile = strcat(directory,'evidence.txt');

%note: this command will only work with matlab R2013 (dataset function not
%included in earlier versions)
msms = dataset2struct(dataset('File',msmsfile));
peptides = dataset2struct(dataset('File',peptidesfile));
proteingroup = dataset2struct(dataset('File',proteinGroupsfile));
evidence = dataset2struct(dataset('File',evidencefile));

[protein] = getProteinStructFromMaxQuantOutput(proteingroup,peptides,msms,evidence,fasta);

% determine if protein is a contaminant
for i = 1:length(protein)
    protein(i).contaminant =  isProteinAContaminant(protein(i));
end



end

