%s2afc_bf1_hier_regress_P300.m - Implements Bayesian hierarchical regression model (Random effects)
%
% Copyright (C) 2020 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  10/23/20        Michael Nunez        Adapted from s2afc_BF1_hier_regress.m

%% Initial

if exist('~/data10') == 7
    data10loc = '~/data10';
elseif exist('/data10') == 7
    data10loc = '/data10';
end
youngloc = [data10loc,'/Stroke2AFC/young'];

saveloc = [data10loc,'/Stroke2AFC/jagsout/'];

load(sprintf('%s/p300peak_svd.mat',youngloc));
load(sprintf('%s/median_RT_DT_thirds_NEWmodel.mat',youngloc));

DTmedian_as = [DTmedian.as.fastest DTmedian.as.medfast DTmedian.as.slowest];
DTmedian_eo = [DTmedian.eo.fastest DTmedian.eo.medfast DTmedian.eo.slowest];
RTmedian_as = [RTmedian.as.fastest RTmedian.as.medfast RTmedian.as.slowest];
RTmedian_eo = [RTmedian.eo.fastest RTmedian.eo.medfast RTmedian.eo.slowest];
p300peaks_AS = [p300peak.as.fastest p300peak.as.medfast p300peak.as.slowest];
p300peaks_EO = [p300peak.eo.fastest p300peak.eo.medfast p300peak.eo.slowest];

participant = [(1:14) (1:14) (1:14)];
P = max(participant);

%% JAGS code for simple linear regression
model = {
'model {'
    '##########'
    '#Prior of Regression intercept (remaining reaction time)'
    'intercept ~ dnorm(.3, pow(.25,-2))'
    ''
    '#Prior of effect of RP latency on reaction time'
    'beta1 ~ dnorm(1,pow(3,-2)) #''Informative'' prior for Bayes Factor calculation'
    ''
    '#Prior of Variability in reaction time'
    'tersubsd ~ dgamma(.3,1)'
    ''
    '#Prior of variance across random intercepts'
    'intersd ~ dgamma(.3, 1)'
    ''
    '#Prior of variance across random effects'
    'beta1sd ~ dgamma(1, 1)'
    ''
    '#Random effects for each participant'
    'for (p in 1:P) {'
        '# Regression intercept'
        'interceptpart[p] ~ dnorm(intercept, pow(intersd,-2))'
        '# Slope intercept'
        'beta1part[p] ~ dnorm(beta1, pow(beta1sd,-2))'
    '}'
    '##########'
    '#Linear regression (normal likelihood)'
    'for (i in 1:N) {'
        'rtpercentile[i] ~ dnorm(interceptpart[participant[i]] + P300peaklat[i]*beta1part[participant[i]], pow(tersubsd,-2) )'
    '}'
'}'
};

%% Code for Trinity

% Track these variables
params = {'intercept', 'beta1', 'tersubsd', 'intersd', 'beta1sd', 'interceptpart', 'beta1part'};

initstruct = @()struct(...
'intercept', rand, ...
'beta1', rand*4 - 2, ...
'tersubsd', rand*.09 + .01,...
'intersd', rand*.19 + .01, ...
'beta1sd', rand*3, ...
'interceptpart', rand(1,P), ...
'beta1part', rand(1,P)*4 - 2);

% Organize data together
dataRTas.rtpercentile = RTmedian_as/1000; %convert to secs
dataRTas.P300peaklat = p300peaks_AS/1000; %convert to secs
dataRTas.N = length(dataRTas.rtpercentile);
dataRTas.participant = participant;
dataRTas.P = P;
dataRTeo.rtpercentile = RTmedian_eo/1000; %convert to secs
dataRTeo.P300peaklat = p300peaks_EO/1000; %convert to secs
dataRTeo.N = length(dataRTeo.rtpercentile);
dataRTeo.participant = participant;
dataRTeo.P = P;
dataDTas.rtpercentile = DTmedian_as/1000; %convert to secs
dataDTas.P300peaklat = p300peaks_AS/1000; %convert to secs
dataDTas.N = length(dataDTas.rtpercentile);
dataDTas.participant = participant;
dataDTas.P = P;
dataDTeo.rtpercentile = DTmedian_eo/1000; %convert to secs
dataDTeo.P300peaklat = p300peaks_EO/1000; %convert to secs
dataDTeo.N = length(dataDTeo.rtpercentile);
dataDTeo.participant = participant;
dataDTeo.P = P;

%% Run P300 peak vs reaction time in action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_hier_regress_P300_RTas';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats_RTas, chains_RTas, diagnostics_RTas, info_RTas] = callbayes('jags', ...
    'model', model, ...
    'data', dataRTas, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_RT] = ksdensity(chains_RTas.beta1(:),xi);
% figure;
% plot(xi,density_RT,'r');
% line(xi,normpdf(xi,1,3));
% title('P300 peak versus RT Percentiles')
numerator = density_RT(xi==1);
denominator = normpdf(1,1,3);
bf_RTas = numerator/denominator;
hierinter_RTas = prctile(chains_RTas.intercept(:),[2.5,50,97.5]);
hierslope_RTas = prctile(chains_RTas.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the P300 vs RT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',hierinter_RTas(2),hierinter_RTas(1),hierinter_RTas(3));
fprintf('The hierarchical slope parameter of the P300 vs RT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',hierslope_RTas(2),hierslope_RTas(1),hierslope_RTas(3));
fprintf('The Bayes Factor of hierarchical slope parameter is %.3f ! \n',bf_RTas);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_RTas', 'chains_RTas', 'diagnostics_RTas','info_RTas','params','dataRTas','model','bf_RTas');
system('rm -r wdir');

%% Run P300 peak vs decision time in action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_hier_regress_P300_DTas';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats_DTas, chains_DTas, diagnostics_DTas, info_DTas] = callbayes('jags', ...
    'model', model, ...
    'data', dataDTas, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_DT] = ksdensity(chains_DTas.beta1(:),xi);
% figure;
% plot(xi,density_DT,'r');
% line(xi,normpdf(xi,1,3));
% title('P300 peak versus RT Percentiles')
numerator = density_DT(xi==1);
denominator = normpdf(1,1,3);
bf_DTas = numerator/denominator;
hierinter_DTas = prctile(chains_DTas.intercept(:),[2.5,50,97.5]);
hierslope_DTas = prctile(chains_DTas.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the P300 vs DT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',hierinter_DTas(2),hierinter_DTas(1),hierinter_DTas(3));
fprintf('The hierarchical slope parameter of the P300 vs DT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',hierslope_DTas(2),hierslope_DTas(1),hierslope_DTas(3));
fprintf('The Bayes Factor of hierarchical slope parameter is %.3f ! \n',bf_DTas);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_DTas', 'chains_DTas', 'diagnostics_DTas','info_DTas','params','dataDTas','model','bf_DTas');
system('rm -r wdir');

%% Run P300 peak vs reaction time in execution only (EO)action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_hier_regress_P300_RTeo';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats_RTeo, chains_RTeo, diagnostics_RTeo, info_RTeo] = callbayes('jags', ...
    'model', model, ...
    'data', dataRTeo, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_RT] = ksdensity(chains_RTeo.beta1(:),xi);
% figure;
% plot(xi,density_RT,'r');
% line(xi,normpdf(xi,1,3));
% title('P300 peak versus RT Percentiles')
numerator = density_RT(xi==1);
denominator = normpdf(1,1,3);
bf_RTeo = numerator/denominator;
hierinter_RTeo = prctile(chains_RTeo.intercept(:),[2.5,50,97.5]);
hierslope_RTeo = prctile(chains_RTeo.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the P300 vs RT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',hierinter_RTeo(2),hierinter_RTeo(1),hierinter_RTeo(3));
fprintf('The hierarchical slope parameter of the P300 vs RT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',hierslope_RTeo(2),hierslope_RTeo(1),hierslope_RTeo(3));
fprintf('The Bayes Factor of hierarchical slope parameter is %.3f ! \n',bf_RTeo);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_RTeo', 'chains_RTeo', 'diagnostics_RTeo','info_RTeo','params','dataRTeo','model','bf_RTeo');
system('rm -r wdir');

%% Run P300 peak vs decision time in execution only (EO) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_hier_regress_P300_DTeo';

modelname = [modeltype,timestr];

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 0; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats_DTeo, chains_DTeo, diagnostics_DTeo, info_DTeo] = callbayes('jags', ...
    'model', model, ...
    'data', dataDTeo, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

%Calculate Bayes Factor using Savage-Dickey ratio
xi = -4:.01:4;
[density_DT] = ksdensity(chains_DTeo.beta1(:),xi);
% figure;
% plot(xi,density_DT,'r');
% line(xi,normpdf(xi,1,3));
% title('P300 peak versus RT Percentiles')
numerator = density_DT(xi==1);
denominator = normpdf(1,1,3);
bf_DTeo = numerator/denominator;
hierinter_DTeo = prctile(chains_DTeo.intercept(:),[2.5,50,97.5]);
hierslope_DTeo = prctile(chains_DTeo.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the P300 vs DT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',hierinter_DTeo(2),hierinter_DTeo(1),hierinter_DTeo(3));
fprintf('The hierarchical slope parameter of the P300 vs DT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',hierslope_DTeo(2),hierslope_DTeo(1),hierslope_DTeo(3));
fprintf('The Bayes Factor of hierarchical slope parameter is %.3f ! \n',bf_DTeo);


save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_DTeo', 'chains_DTeo', 'diagnostics_DTeo','info_DTeo','params','dataDTeo','model','bf_DTeo');
system('rm -r wdir');