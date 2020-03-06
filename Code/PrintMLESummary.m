%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DisplayMLESummary.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% * Display summary of MLE.

function PrintMLESummary(mleResultStruct, filepath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected the following inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mleResultStruct: Structure generated by MLEEstimation().
% filepath: String path to output file.

if nargin < 2
   error('Require mleResultStruct, filepath args.') 
end

paramNames = mleResultStruct.paramNames;
estimators = mleResultStruct.estimators;
likelihood = mleResultStruct.likelihood;
nulls = mleResultStruct.nulls;
pvals = mleResultStruct.pvals; 
tvals = mleResultStruct.t_obs;
stderrs = mleResultStruct.stderrs;
hypStruct = mleResultStruct.hypothesisparams;

file = fopen(filepath, 'w');
fprintf(file, "-----------------------\n");
fprintf(file, "MLE Estimation Results:\n");
fprintf(file, "-----------------------\n");
fprintf(file, strcat("Log-Likelihood:,",num2str(likelihood),"\n"));
fprintf(file, strcat("alpha:,", num2str(hypStruct.alpha),"\n"));

% Print headers:
fprintf(file, "Param:,Value:,Hypothesis:,T-Obs:,StdErr:,PVal:\n");
% Print data:
formatLine = "%s,%f,%f,%f,%f,%f\n";
for est = 1:length(estimators)
    fprintf(file, formatLine, paramNames(est), estimators(est), nulls(est),tvals(est),stderrs(est),pvals(est));
end

fclose(file);

end