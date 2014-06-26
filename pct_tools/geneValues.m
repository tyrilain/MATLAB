% GENE VALUES AggreateS the values of a gene over several time cohorts into an array.
% 
%    VALUES = GENEVALUES(VPC, GENENAME, COHORTINDICES) returns an array that has all
%    the expression values of GENENAME in the virtual point cloud VPC for the
%    cohorts with indices COHORTINDICES
%    
%    If no values are available for a particular cohort NaN's are output
%    instead
%
% Zeba Wunderlich, 2009-08-17

function values = geneValues(vpc, geneName, cohortIndices)


values = [];
for k=1:length(cohortIndices)
    fieldString = [geneName '__' num2str(cohortIndices(k))];
    if (~isfield(vpc, fieldString))
        values(end+1:end+vpc.metadata.nuclear_count) = 0;
    else
        eval(['values(end+1:end+vpc.metadata.nuclear_count) = vpc.' fieldString ';']);
    end
end
values = values';

    