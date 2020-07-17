function trials = organizetrials(artifact,blockorder,correctresponses,responseside,responsetime);
%organizes list of trials for use in analysis 
%INPUT: artifact - artifact matrix from cleaned data (comes as 1 is bad)
%	blockorder - indexes which trial type
%	correctresponses - indexes which stimulus, 1 = left, 2 = right
%	responseside - subject response, 1 = left, 2 = right
%	responsetime - reaction time
if nargin < 5
exit('insufficient input arguments')
end;
badtrials = find(sum(artifact) == size(artifact,1));
badresponsetrials = find(responsetime > 2399);
noresponsetrials = find(isnan(responsetime));
trials.badtrials = badtrials;
trials.badresponsetrials = badresponsetrials;
trials.noresponsetrials = noresponsetrials;
badtrials = unique([badtrials badresponsetrials noresponsetrials]);
eo.left = setdiff(find(blockorder == 1),badtrials);
eo.right = setdiff(find(blockorder == 2),badtrials);
as.left = setdiff(find(blockorder == 3 & responseside == 1),badtrials);
as.right = setdiff(find(blockorder == 3 & responseside == 2),badtrials);

testaccuracy = (responseside == correctresponses);

correct.eo.left = eo.left(find(testaccuracy(eo.left)));
correct.eo.right = eo.right(find(testaccuracy(eo.right)));
correct.as.left = as.left(find(testaccuracy(as.left)));
correct.as.right = as.right(find(testaccuracy(as.right)));
%organize output
trials.eo = eo;
trials.as = as;
trials.correct = correct;