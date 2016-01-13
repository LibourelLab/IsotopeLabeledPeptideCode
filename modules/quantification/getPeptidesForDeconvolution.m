function [peps_for_deconv] = getPeptidesForDeconvolution(peptides)


peps_for_deconv = peptides(cellfun(@length,{peptides.feature}) == 1);

% for each peptide, append field with appropriate field names
for i = 1:length(peps_for_deconv)
    peps_for_deconv(i).measured = peps_for_deconv(i).feature.ria;
    peps_for_deconv(i).isotope =  peps_for_deconv(i).feature.isotope;
    peps_for_deconv(i).weight =  sum([peps_for_deconv(i).feature.ions(:)]);
end

end