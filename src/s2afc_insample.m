%S2AFC_INSAMPLE - Finds posterior predictive distributions of model fits in s2afc_behavmodel.m

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  10/14/20        Michael Nunez               Original code (pdmfinal_inout_regression_N200.py reference)
%  10/20/20        Michael Nunez           Include scatter plots and regression lines
%  10/22/20        Michael Nunez      Calculate pearson Rho between true 10 percentiles and NDT


%% Initial
if exist('~/data10') == 7
    data10loc = '~/data10';
elseif exist('/data10') == 7
    data10loc = '/data10';
end

jagsoutloc = [data10loc,'/Stroke2AFC/jagsout/jagsmodelbasic_lapse_8_17_16_51.mat'];

jagsoutall = load(jagsoutloc);
jagsout = readjagsout(jagsoutall.stats,jagsoutall.diagnostics);
jagsin = jagsoutall.data;

%% Generate posterior predictive distributions

%Number of chains
nchains = size(jagsoutall.chains.alphasub_1_1,2);

%Number of samples per chain
nsamps = size(jagsoutall.chains.alphasub_1_1,1);

%Number of samples in total
ntotalsamps = nchains*nsamps;

%Set seed
rng(2020); 

N = jagsin.N;
Niter = 100;

fprintf('Generating %d * %d = %d samples...\n',N,Niter,N*Niter);
correctRT_pred = zeros(N, Niter);
Alpha = zeros(N, Niter);
Ndt = zeros(N, Niter);
Delta = zeros(N, Niter);
problapse = zeros(N, Niter);
for n=1:N,
    randsamps = randi(ntotalsamps,Niter);
    whichchain = ceil(randsamps/nsamps);
    whichsamp = mod(randsamps,nsamps)+1; %Note that this indexing does not exactly match to an ordering of "randsamps" but that doesn't matter as long as the samples are unique
    condindex = jagsin.condition(n); %1 - AS, 2 - EO
    subindex = jagsin.subject(n);
    for i=1:Niter,
        %Use known parameters for each participant to generate choice-RTs
        Alpha(n,i) = jagsoutall.chains.(sprintf('alphasub_%d_%d',condindex,subindex))(whichsamp(i),whichchain(i));
        Ndt(n,i) = jagsoutall.chains.(sprintf('tersub_%d_%d',condindex,subindex))(whichsamp(i),whichchain(i));
        Delta(n,i) = jagsoutall.chains.(sprintf('deltasub_%d_%d',condindex,subindex))(whichsamp(i),whichchain(i));
        %Find lapse probability for this participant
        problapse(n,i) = jagsoutall.chains.(sprintf('probsub_%d_%d_2',condindex,subindex))(whichsamp(i),whichchain(i));

        notlapse = (problapse(n,i) < rand);
        if notlapse,
            [rt, acc] = simuldiff([Alpha(n,i)/10, Ndt(n,i), 0, (Alpha(n,i)/10)*.5, 0, 0, Delta(n,i)/10],1); %Note simuldiff uses absolute start point and diffusion coefficient =.1
            correctRT_pred(n,i) = (acc*2 - 1)*rt; % negative RTs for incorrect responses
        else
            correctRT_pred(n,i) = rand*jagsin.maxrt(condindex,subindex)*2 - jagsin.maxrt(condindex,subindex); %range of (-maxrt, maxrt)
        end
    end
end

%% Evaluate in-sample prediction
% See: Ratcliff, R., Smith, P. L., Brown, S. D., & McKoon, G. (2016). Diffusion decision model: Current issues and history. Trends in cognitive sciences, 20(4), 260-281.

nplotconds = jagsin.nconds*2; % Double the number of conditions for errors and correct responses
percentiles = [10, 30, 50, 70, 90];
rtpercentiles_true = nan(length(percentiles),jagsin.nsubs,nplotconds);
rtpercentiles_pred = nan(length(percentiles),jagsin.nsubs,nplotconds);
accuracy_true = nan(jagsin.nsubs,jagsin.nconds);
accuracy_pred = nan(jagsin.nsubs,jagsin.nconds);


for c=1:jagsin.nconds,
    for s=1:jagsin.nsubs,
        whereboth = ((jagsin.subject == s) & (jagsin.condition == c));
        wherecorrect = (jagsin.y > 0) & whereboth;
        whereerror = (jagsin.y < 0) & whereboth;
        accuracy_true(s,c) = mean((jagsin.y(whereboth) > 0)); %True accuracies
        accuracy_pred(s,c) = mean(mean((correctRT_pred(whereboth,:) > 0),1),2); %Predicted accuracies
        for r=1:length(percentiles),
            rtpercentiles_true(r,s,c) = prctile(jagsin.y(wherecorrect),percentiles(r)); %Correct RTs
            rtpercentiles_true(r,s,c+2) = prctile(abs(jagsin.y(whereerror)),percentiles(r)); %Error RTs
            all_predvals = reshape(correctRT_pred(whereboth,:),1,sum(whereboth)*Niter); %All relevant predicted values collapsed to a vector
            rtpercentiles_pred(r,s,c) = prctile(all_predvals(all_predvals > 0),percentiles(r)); %Correct predicted RTs
            rtpercentiles_pred(r,s,c+2) = prctile(abs(all_predvals(all_predvals < 0)),percentiles(r)); %Incorrect predicted RTs
        end
    end
end

%Print mean accuracies and ranges
fprintf('The mean accuracies of participants in the AS task was %.3f with a min of %.3f and max of %.3f! \n',mean(accuracy_true(:,1)),min(accuracy_true(:,1)),max(accuracy_true(:,1)));
fprintf('The mean accuracies of participants in the EO task was %.3f with a min of %.3f and max of %.3f! \n',mean(accuracy_true(:,2)),min(accuracy_true(:,2)),max(accuracy_true(:,2)));


%Compute R^2 of prediction values
accuracy_rsquared = nan(jagsin.nconds,1);
rtpercentiles_rsquared = nan(length(percentiles),nplotconds);

for c=1:jagsin.nconds,
    divisor = jagsin.nsubs -1;
    %Mean squared error of prediction
    MSEP_acc = sum((accuracy_true(:,c) - accuracy_pred(:,c)).^2) / divisor;
    %Variance estimate of the true values
    vartrue_acc = sum((accuracy_true(:,c) - mean(accuracy_true(:,c)) ).^2) / divisor;
    accuracy_rsquared(c) = 1 - (MSEP_acc / vartrue_acc);
    for r=1:length(percentiles),
        divisor = sum(isfinite(squeeze(rtpercentiles_true(r,:,c)))) - 1;
        %Mean squared error of prediction
        MSEP_rt = nansum((squeeze(rtpercentiles_true(r,:,c)) - squeeze(rtpercentiles_pred(r,:,c)) ).^2) / divisor;
        %Variance estimate of the true values
        vartrue_rt = nansum((squeeze(rtpercentiles_true(r,:,c)) - nanmean(squeeze(rtpercentiles_true(r,:,c))) ).^2) / divisor;
        rtpercentiles_rsquared(r,c) = 1 - (MSEP_rt / vartrue_rt);
        %Mean squared error of prediction
        MSEP_rt = nansum((squeeze(rtpercentiles_true(r,:,c+2)) - squeeze(rtpercentiles_pred(r,:,c+2)) ).^2) / divisor;
        %Variance estimate of the true values
        vartrue_rt = nansum((squeeze(rtpercentiles_true(r,:,c+2)) - nanmean(squeeze(rtpercentiles_true(r,:,c+2))) ).^2) / divisor;
        rtpercentiles_rsquared(r,c+2) = 1 - (MSEP_rt / vartrue_rt);
    end
end

%% Make table of R^2_pred values
addpath('../external'); %Find latextable.m
Vertlabels = {'\textbf{10th RT Percentile}', '\textbf{30th RT Percentile}', '\textbf{50th RT Percentile}', '\textbf{70th RT Percentile}', '\textbf{90th RT Percentile}'};
Horizlabels = {'\textbf{Action Selection (AS)}','\textbf{Execution Only (EO)}'};
tableinput = {sprintf('%.0f \\%%', rtpercentiles_rsquared(1, 1)*100), sprintf('%.0f \\%%', rtpercentiles_rsquared(1, 2)*100); ...
            sprintf('%.0f \\%%', rtpercentiles_rsquared(2, 1)*100), sprintf('%.0f \\%%', rtpercentiles_rsquared(2, 2)*100); ...
            sprintf('%.0f \\%%', rtpercentiles_rsquared(3, 1)*100), sprintf('%.0f \\%%', rtpercentiles_rsquared(3, 2)*100); ...
            sprintf('%.0f \\%%', rtpercentiles_rsquared(4, 1)*100), sprintf('%.0f \\%%', rtpercentiles_rsquared(4, 2)*100); ...
            sprintf('%.0f \\%%', rtpercentiles_rsquared(5, 1)*100), sprintf('%.0f \\%%', rtpercentiles_rsquared(5, 2)*100)};


filesave = '../../Paper/CoBB/Supplementals/Insample_table.tex';
latextable(tableinput, 'Horiz', Horizlabels, 'Vert', Vertlabels, 'Hline', [0,1,NaN],'name',filesave);


%% Calculate pearson correlation between true 10 percentiles and NDT in each task

Ndt = nan(jagsin.nconds,jagsin.nsubs);

for c=1:jagsin.nconds,
    for s=1:jagsin.nsubs,
        %Find medians of subject-level parameters
        Ndt(c,s) = prctile(jagsoutall.chains.(sprintf('tersub_%d_%d',c,s))(:),50);
    end
end

[rhoAS, pvalAS] = corr(squeeze(rtpercentiles_true(1,:,1))',Ndt(1,:)');
fprintf('The Pearson correlation coefficient between 10th RT percentiles and NDT is %.3f (p = %.3f) in the AS task! \n',rhoAS,pvalAS);
[rhoEO, pvalEO] = corr(squeeze(rtpercentiles_true(1,:,2))',Ndt(2,:)');
fprintf('The Pearson correlation coefficient between 10th RT percentiles and NDT is %.3f (p = %.3f) in the EO task! \n',rhoEO,pvalEO);