%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ConvertSymbolicFunction.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert symbolic function with parameter
% set x_i (data), a-z (parameters) for use
% with MLEEstimation functions.

function func_handle = ConvertSymbolicFunction(symbolFunc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expecting following parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbolStr = string(symbolFunc);
% Convert all x_i parameters, abcd etc parameters, to element of vector:
matches = unique(regexp(symbolStr, "x_[0-9]+", "match"));
for i = 1:length(matches)
    match = matches(i);
    val = regexp(match, "[0-9]+", "match");
    symbolStr = regexprep(symbolStr, match, strcat('x(', val, ')'));
end
matches = unique(regexp(symbolStr, "p_[0-9]+", "match"));
for i = 1:length(matches)
    match = matches(i);
    val = regexp(match, "[0-9]+", "match");
    symbolStr = regexprep(symbolStr, match, strcat('params(', val, ')'));
end

symbolStr = strcat('@(x, params)', symbolStr);
func_handle = str2func(symbolStr);

end

function mapped = MapChar(charStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expected inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charStr: string containing characters denoting parameter names in
% symbolic function string. Assumes that parameters are denoted by
% a-z. 

outVal = double(char(charStr)) - double('a') + 1;
mapped = num2str(outVal);

end