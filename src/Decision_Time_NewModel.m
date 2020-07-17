%% Calculates the decision-time based on the jags outputs this is done for each participant
clear all 
load('/home/kitty/data10/Stroke2AFC/jagsout/jagsmodelbasic_lapse_8_17_16_51.mat')
jagsout = readjagsout(stats,diagnostics);
nondt.as = jagsout.median(57:2:84);
nondt.eo = jagsout.median(58:2:84);
nondt.all = jagsout.median(57:84);
nondt.as = nondt.as*1000;
nondt.eo = nondt.eo*1000;

data = load('subjid_cleaned.mat');

%find all the good trials w. goodbad function

[goodchans,goodtrials,badchans,badtrials] = goodbad(data);

% get all good trials (excludes bad trials & those with NaN in RTs)
badresponsetrials = find(data.responsetime > 2399 | data.responsetime < 200);
noresponsetrials = find(isnan(data.responsetime));
longtrials = find(data.responsetime > 2399);
usabletrials = setdiff(goodtrials,noresponsetrials);
usabletrials = setdiff(usabletrials,badresponsetrials);

as_trials = intersect(find(data.blockorder == 3), usabletrials);
eo_trials = setdiff(usabletrials,as_trials);

as_dt = data.responsetime(as_trials) - nondt.as(14); %make sure to change the nondt.as(15) index according to what number the subject is
eo_dt = data.responsetime(eo_trials) - nondt.eo(14);

as_trialsDT = [as_trials(:),as_dt(:)];
eo_trialsDT = [eo_trials(:),eo_dt(:)];

%% everything after this is for splitting the DT into tertiles
as_cutoff = prctile(as_dt,[33 66]);
eo_cutoff = prctile(eo_dt,[33 66]);

for j = 1:length(as_trials)
    if as_trialsDT(j,2) < as_cutoff(1)
        trials.as.fastest(j) = as_trialsDT(j);
    elseif as_trialsDT(j,2) > as_cutoff(2)
         trials.as.slowest(j) = as_trialsDT(j);
    else trials.as.medfast(j) = as_trialsDT(j);
    end
end

trials.as.fastest(trials.as.fastest == 0) = [ ]; 
trials.as.medfast(trials.as.medfast == 0) = [ ];
trials.as.slowest(trials.as.slowest == 0) = [ ]; 

for j = 1:length(eo_trials)
    if eo_trialsDT(j,2) < eo_cutoff(1)
        trials.eo.fastest(j) = eo_trialsDT(j);
    elseif eo_trialsDT(j,2) > eo_cutoff(2)
         trials.eo.slowest(j) = eo_trialsDT(j);
    else trials.eo.medfast(j) = eo_trialsDT(j);
    end
end

trials.eo.fastest(trials.eo.fastest == 0) = [ ]; 
trials.eo.medfast(trials.eo.medfast == 0) = [ ];
trials.eo.slowest(trials.eo.slowest == 0) = [ ]; 

findfast_as = ismember(as_trialsDT(:,1),trials.as.fastest);
% dt.as.fastest = find(findfast_as);
dt.as.fastest = as_trialsDT(findfast_as,2);

findmed_as = ismember(as_trialsDT(:,1),trials.as.medfast);
% dt.as.medfast = find(findmed_as);
dt.as.medfast = as_trialsDT(findmed_as,2);

findslow_as = ismember(as_trialsDT(:,1),trials.as.slowest);
% dt.as.slowest = find(findslow_as);
dt.as.slowest = as_trialsDT(findslow_as,2);

findfast_eo = ismember(eo_trialsDT(:,1),trials.eo.fastest);
% dt.eo.fastest = find(findfast_eo);
dt.eo.fastest = eo_trialsDT(findfast_eo,2);

findmed_eo = ismember(eo_trialsDT(:,1),trials.eo.medfast);
% dt.eo.medfast = find(findmed_eo);
dt.eo.medfast = eo_trialsDT(findmed_eo,2);

findslow_eo = ismember(eo_trialsDT(:,1),trials.eo.slowest);
% dt.eo.slowest = find(findslow_eo);
dt.eo.slowest = eo_trialsDT(findslow_eo,2);

save subjid_dt_thirds_NEW.mat trials dt as_trialsDT eo_trialsDT