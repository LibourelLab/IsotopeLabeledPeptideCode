function [ mass_traces ] = getMassTraces(mass_trace_text_file,rawdata)
%getMassTraces is a function imports the text file generated from the
%python script, reduces the outut and re-maps the data back to raw data

%input mass traces into a data structure
mass_traces =feature2struct(mass_trace_text_file);
%map features to raw data 
mass_traces = mapFeatureData2RawData(rawdata,mass_traces,true);
%reduce mapped features to get mean retention time and 
mass_traces = reduceMappedFeatures(mass_traces);


end

