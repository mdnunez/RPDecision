%% this script splits up the RT into tertiles for each subject and saves out a mat file
%Figure 5b, Figure 6 (RT split averages)
clear all
%load subject data in:

data = load('subjid_cleaned.mat');

%find all the good trials w. goodbad function

[goodchans,goodtrials,badchans,badtrials] = goodbad(data);

% get all good trials (excludes bad trials & those with NaN in RTs)
badresponsetrials = find(data.responsetime > 1600 | data.responsetime < 200);
noresponsetrials = find(isnan(data.responsetime));
longtrials = find(data.responsetime > 1600);
usabletrials = setdiff(goodtrials,noresponsetrials);
usabletrials = setdiff(usabletrials,badresponsetrials);

as_trials = intersect(find(data.blockorder == 3), usabletrials);
eo_trials = setdiff(usabletrials,as_trials);

as_responsetimes = data.responsetime(as_trials);
eo_responsetimes = data.responsetime(eo_trials);

as_trialsRT = [as_trials(:),as_responsetimes(:)];
eo_trialsRT = [eo_trials(:),eo_responsetimes(:)];

as_cutoff = prctile(as_responsetimes,[33 66]);
eo_cutoff = prctile(eo_responsetimes,[33 66]);

for j = 1:length(as_trials)
    if as_trialsRT(j,2) < as_cutoff(1)
        trials.as.fastest(j) = as_trialsRT(j);
    elseif as_trialsRT(j,2) > as_cutoff(2)
            trials.as.slowest(j) = as_trialsRT(j);
    else trials.as.medfast(j) = as_trialsRT(j);
    end
end

trials.as.fastest(trials.as.fastest == 0) = [ ]; 
trials.as.medfast(trials.as.medfast == 0) = [ ];
trials.as.slowest(trials.as.slowest == 0) = [ ]; 

for j = 1:length(eo_trials)
    if eo_trialsRT(j,2) < eo_cutoff(1)
        trials.eo.fastest(j) = eo_trialsRT(j);
    elseif eo_trialsRT(j,2) > eo_cutoff(2)
            trials.eo.slowest(j) = eo_trialsRT(j);
    else trials.eo.medfast(j) = eo_trialsRT(j);
    end
end

trials.eo.fastest(trials.eo.fastest == 0) = [ ]; 
trials.eo.medfast(trials.eo.medfast == 0) = [ ];
trials.eo.slowest(trials.eo.slowest == 0) = [ ]; 




save subjid_RTthirds.mat trials 