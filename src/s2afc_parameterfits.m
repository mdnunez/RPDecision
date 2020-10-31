%S2AFC_PARAMETERFITS - Generates table of parameter fits of model in s2afc_behavmodel.m

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  10/19/20        Michael Nunez               Original code



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

%% Find percentiles of posterior distributions

N = jagsin.N;

percentiles = [2.5, 25, 50, 75, 97.5];

fprintf('Finding medians and 95%% credible intervals for each parameter...\n');
Alpha = nan(length(percentiles),jagsin.nconds,jagsin.nsubs);
Ndt = nan(length(percentiles),jagsin.nconds,jagsin.nsubs);
Delta = nan(length(percentiles),jagsin.nconds,jagsin.nsubs);
problapse = nan(length(percentiles),jagsin.nconds,jagsin.nsubs);

Alphacond = nan(length(percentiles),jagsin.nconds);
Ndtcond = nan(length(percentiles),jagsin.nconds);
Deltacond = nan(length(percentiles),jagsin.nconds);

for c=1:jagsin.nconds,
    %Find percentiles of hierarchical condition-level parameters
    Alphacond(:,c) = prctile(jagsoutall.chains.(sprintf('alphacond_%d',c))(:),percentiles);
    Ndtcond(:,c) = prctile(jagsoutall.chains.(sprintf('tercond_%d',c))(:),percentiles);
    Deltacond(:,c) = prctile(jagsoutall.chains.(sprintf('deltacond_%d',c))(:),percentiles);
    for s=1:jagsin.nsubs,
        %Find percentiles of subject-level parameters
        Alpha(:,c,s) = prctile(jagsoutall.chains.(sprintf('alphasub_%d_%d',c,s))(:),percentiles);
        Ndt(:,c,s) = prctile(jagsoutall.chains.(sprintf('tersub_%d_%d',c,s))(:),percentiles);
        Delta(:,c,s) = prctile(jagsoutall.chains.(sprintf('deltasub_%d_%d',c,s))(:),percentiles);
        %Find lapse probability for this participant
        problapse(:,c,s) = prctile(jagsoutall.chains.(sprintf('probsub_%d_%d_2',c,s))(:),percentiles);
    end
end

%% Generate LaTeX table
addpath('../external'); %Find latextable.m
Vertlabels = {'\textbf{$\tau$ (sec) }','\textbf{$\alpha$ (evid.)}', '\textbf{$\delta$ (evid./sec)}', '\textbf{$\lambda$ (probability)}'};
Horizlabels = {'\textbf{AS Mean (95\% CI)}','\textbf{AS Part. Range}', '\textbf{EO Mean (95\% CI)}', '\textbf{EO Part. Range}'};
tableinput = {sprintf('%.2f [%.2f, %.2f]', Ndtcond(3,1), Ndtcond(1,1), Ndtcond(5,1)),...
                sprintf('%.2f to %.2f', min(squeeze(Ndt(3,1,:))), max(squeeze(Ndt(3,1,:))) ),...
                sprintf('%.2f [%.2f, %.2f]', Ndtcond(3,2), Ndtcond(1,2), Ndtcond(5,2)),...
                sprintf('%.2f to %.2f', min(squeeze(Ndt(3,2,:))), max(squeeze(Ndt(3,2,:))) ); ...
                sprintf('%.2f [%.2f, %.2f]', Alphacond(3,1), Alphacond(1,1), Alphacond(5,1)),...
                sprintf('%.2f to %.2f', min(squeeze(Alpha(3,1,:))), max(squeeze(Alpha(3,1,:))) ),...
                sprintf('%.2f [%.2f, %.2f]', Alphacond(3,2), Alphacond(1,2), Alphacond(5,2)),...
                sprintf('%.2f to %.2f', min(squeeze(Alpha(3,2,:))), max(squeeze(Alpha(3,2,:))) ); ...
                sprintf('%.2f [%.2f, %.2f]', Deltacond(3,1), Deltacond(1,1), Deltacond(5,1)),...
                sprintf('%.2f to %.2f', min(squeeze(Delta(3,1,:))), max(squeeze(Delta(3,1,:))) ),...
                sprintf('%.2f [%.2f, %.2f]', Deltacond(3,2), Deltacond(1,2), Deltacond(5,2)),...
                sprintf('%.2f to %.2f', min(squeeze(Delta(3,2,:))), max(squeeze(Delta(3,2,:))) ); ...
                sprintf('%.2f\\%% **', mean(squeeze(problapse(3,1,:)))*100 ), ...
                sprintf('%.2f\\%% to %.2f\\%%', min(squeeze(problapse(3,1,:)))*100, max(squeeze(problapse(3,1,:)))*100 ),...
                sprintf('%.2f\\%% **', mean(squeeze(problapse(3,2,:)))*100 ), ...
                sprintf('%.2f\\%% to %.2f\\%%', min(squeeze(problapse(3,2,:)))*100, max(squeeze(problapse(3,2,:)))*100 ) };

filesave = '../../Paper/CoBB/Supplementals/Parameter_table.tex';
latextable(tableinput, 'Horiz', Horizlabels, 'Vert', Vertlabels, 'Hline', [0,1,NaN],'name',filesave);