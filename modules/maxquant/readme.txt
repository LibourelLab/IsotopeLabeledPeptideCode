%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxQuant Data Parser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Created 6-19-2013;
Edited 6-19-2013;

Requirements 
1) Directory that contains all text output from maxquant searches
2) Fasta file from unlabaled experiment
3) Raw data in mzXML format

Functionality includes:
1) File i/o of maxquant search output in matlab R2013 (must include dataset function)
2) fasta file reduction from proteins identified in experiment. Currently only uniprot format is supported.
3) Mapping precursor ion information from measured peptides to peptide data structure using mzxml raw data

*****Suggested Workflow********
1) run maxquant
2) generate protein structure using maxQuantOutput2protein.m
3) reduce fasta file using reduceFasta.m
4) re-run maxquant with new fasta file
5) repeat step 2
6) run mapPrecursorMZtoProteinStruct.m to append corrected m/z value to each peptide

