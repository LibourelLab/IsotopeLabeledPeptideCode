function [ fasta_new ] = reduceFasta(proteinGroups,fastafile,format)
%reduce fasta file with by IDed protein groups


% reduce protein groups to gid rid of contaminants
%proteinGroups=proteinGroups(cellfun(@isempty,{proteinGroups(:).Contaminant},'uniform',1));

fastaString = cellfun(@(x) regexp(x,'\|','split'),{proteinGroups(:).Fastaheaders},'uniform',0);

%6-19-2013: currently only expression is for UNIPROT databases, need to extend for
%NCBI and others...
if (format == 'UniProt')
        expression = '[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]|[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]';
else
end
an = regexp({proteinGroups(:).Fastaheaders},expression,'match');
an = unique([an{:}]);

mapped_index = cellfun(@(x) find(not(cellfun(@isempty,regexp({fastafile(:).Header},x,'match'),'uniform',1))),an,'uniform',0);
mapped_index = cell2mat(mapped_index);
fasta_new  = fastafile(mapped_index);
%%need to add unique statement
end

