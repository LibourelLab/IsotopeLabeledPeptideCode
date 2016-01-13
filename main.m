%% MATLAB Script for running automated pipeline identification and quantifiaction of isotope labeled peptides
% Joshua E. Goldford
% 1-9-2016
% 
% note that the example data files are not provided via GitHub and can be found in the supplement of the following manuscript
% Goldford JE, Libourel IG.  Unsupervised Identification of Isotope Labeled Pepides
%
%
% please use this script along with the document tutorial.pdf

%% ***********************************************************************************************************************
%  ***********************************  PROTEIN DATABASE FROM MAXQUANT  **************************************************

% set currrent directoy to variable
dir = pwd;

% add paths where functions are
addpath('modules/maxquant/lib');
addpath('modules/openms/lib/matlab');
addpath('modules/openms/lib/matlab');
addpath('modules/extraction/lib');
addpath('modules/quantification/lib');
addpath('modules/msfilereader_matlab_api/64/mex');
addpath('modules/msfilereader_matlab_api/64/sub');
addpath('modules/msfilereader_matlab_api/64');


% parse maxqquant output into a matlab structure
protein = maxQuantOutput2protein([dir,'/data/ms/unlabeled/maxquant/'],[dir,'/data/sequence/E_coli_04022013.fasta']);

% map m/z precursor for each identified peptide from maxquant - 
protein = mapPrecursorMZtoProteinStruct(protein,[dir,'/data/ms/unlabeled/mzxml/unlabeled.mzXML']);

% remove contaminats from protein database
protein([protein.contaminant]) = [];



%% ***********************************************************************************************************************
%  ***********************************  IMPORT LABELED DATA  *************************************************************

% Note: this step HAS to be done using the msfilereader-matlab-api.  Please
% look at the documentation on the GitHub.  The import of .RAW data can only
% be done using a Windows machine with MSFileReader installed, as well as
% MATLAB and a c++ compiler.

path_to_raw_data_file = '';
raw = getRawData(path_to_raw_data_file);


%% ***********************************************************************************************************************
%  ***********************************  IMPORT MASS TRACES FROM OPENMS *****************************************


% get mass traces from Opem MS output (see open MS module documentation to
% get a file).  This funciton matches the mass traces to the data imported
% from the raw data function... This is important because the raw data
% contains a lot of useful quantitative information that mzxml files (and
% subsequently all data from OpenMS) do not have.  
mass_traces = getMassTraces([dir,'/data/ms/labeled/mass_traces.tsv'],raw);
clear raw;



%% ***********************************************************************************************************************
%  ***********************************  FILTER MASS TRACES  *****************************************


% filter mass traces (we will just do one protein in this example);  first
% filtering step is for a MAC of 5 ppm and RTW of 1 min
pro_idx = 73;
mtf.mac5rtw1 = filterMassTraceForProtein(mass_traces,protein(pro_idx),5,1);
mtf.mac10rtw2 = filterMassTraceForProtein(mass_traces,protein(pro_idx),10,2);


%% ***********************************************************************************************************************
%  ***********************************  CLUSTER MASS TRACES  *****************************************

peptides.mac5rtw1= cluster_mass_traces(protein(pro_idx),mtf.mac5rtw1{1});
peptides.mac10rtw2= cluster_mass_traces(protein(pro_idx),mtf.mac10rtw2{1});

%% ***********************************************************************************************************************
%  *********************************** ITERATIVE DECONVOLUTION  *****************************************

% (1) MAC 5 and RTW 1: INITIAL AMINO ACID MID ESTIMATES FROM STRICT FILTERING

% keep only peptides with exactly one feature
deconvPeptideSet = getPeptidesForDeconvolution(peptides.mac5rtw1);

% perform deconvolution of unique features to get an amino acid MID
% estimate from the initial mass trace cutoffs
amino_acid_mid = automated_deconv(deconvPeptideSet,1,[]);

% (1) MAC 10 and RTW 2: SIMULATE PMD & REDUCE CANDIDATE FEATURES BASED ON
% PERCENT AGREEMENT BETWEEN MEASUREMENT AND SIMULATED PMD (CUTOFF = 0.5).
% UPDATE AA MID ESTIMATE

% simulate peptide MIDs from larger retention time window
peptides.mac10rtw2 = simulate_peptide_MIDS(peptides.mac10rtw2,amino_acid_mids);

% reduce the number of peptide candidate features based on a percent agreement cutoff of 0.5
peptides_temp = reduce_features(peptides.mac10rtw2,0.5);

% keep only peptides with exactly one feature for devconvolution
deconvPeptideSet = getPeptidesForDeconvolution(peptides_temp);

% perform deconvolution of unique features to get an amino acid MID
% estimate from the initial mass trace cutoffs
amino_acid_mid = automated_deconv(deconvPeptideSet,1,amino_acid_mid);


% (2) REPEAT (1) WITH CUTOFF = 0.75


% update feature simulated peptide MIDS and percent agreement
peptides.mac10rtw2 = simulate_peptide_MIDS(peptides.mac10rtw2,amino_acid_mids);

% reduce the number of peptide candidate features based on a percent agreement cutoff of 0.75
peptides_temp = reduce_features(peptides.mac10rtw2,0.75);

% keep only peptides with exactly one feature for devconvolution
deconvPeptideSet = getPeptidesForDeconvolution(peptides_temp);

% perform deconvolution of unique features to get an amino acid MID
% estimate from the initial mass trace cutoffs
amino_acid_mid = automated_deconv(deconvPeptideSet,1,amino_acid_mid);


% (3) REPEAT (2) WITH CUTOFF = 0.9

% update feature simulated peptide MIDS and percent agreement
peptides.mac10rtw2 = simulate_peptide_MIDS(peptides.mac10rtw2,amino_acid_mids);

% reduce the number of peptide candidate features based on a percent
% agreement cutoff of 0.90
peptides_temp = reduce_features(peptides.mac10rtw2,0.9);

% keep only peptides with exactly one feature for devconvolution
deconvPeptideSet = getPeptidesForDeconvolution(peptides_temp);

% perform deconvolution of unique features to get an amino acid MID
% estimate from the initial mass trace cutoffs
amino_acid_mid = automated_deconv(deconvPeptideSet,1,amino_acid_mid);

% update feature simulated peptide MIDS and percent agreement
peptides.mac10rtw2 = simulate_peptide_MIDS(peptides.mac10rtw2,amino_acid_mids);


% (4) UPDATE SIMULATED PMD FROM BEST AA MID ESTIMATE AND KEEP ONLY PEPTIDES
% WITH ONE CANDIDATE FEATURE

% reduce the number of peptide candidate features based on a percent
% agreement cutoff of 0.90 again
peptides_temp = reduce_features(peptides.mac10rtw2,0.9);

% keep only peptides with exactly one feature for final assignment
final_peptide_set = getPeptidesForDeconvolution(peptides_temp);







