%s2afc_bf1regression_rtthirds_expanded_priors.m - Calculates Bayes Factors for simple regression models with RP durations with less restrictive priors
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
%  10/27/18        Michael Nunez             Adapted from s2afc_BF1regression_rtthirds.m

%% Initial

if exist('~/data10') == 7
    data10loc = '~/data10';
elseif exist('/data10') == 7
    data10loc = '/data10';
end
youngloc = [data10loc,'/Stroke2AFC/young'];

saveloc = [data10loc,'/Stroke2AFC/jagsout/'];

load(sprintf('%s/rp_onsets.mat',youngloc));
% load(sprintf('%s/all_onsets_100.mat',youngloc));
load(sprintf('%s/median_RT_DT_thirds_NEWmodel.mat',youngloc));

DTmedian_as = [DTmedian.as.fastest DTmedian.as.medfast DTmedian.as.slowest];
DTmedian_eo = [DTmedian.eo.fastest DTmedian.eo.medfast DTmedian.eo.slowest];
RTmedian_as = [RTmedian.as.fastest RTmedian.as.medfast RTmedian.as.slowest];
RTmedian_eo = [RTmedian.eo.fastest RTmedian.eo.medfast RTmedian.eo.slowest];
rponsets_as = [abs(rponset.as.fastest) abs(rponset.as.medfast) abs(rponset.as.slowest)];
rponsets_eo = [abs(rponset.eo.fastest) abs(rponset.eo.medfast) abs(rponset.eo.slowest)];
% rponsets_as = [abs(responsets_100.as.fastest) abs(responsets_100.as.medfast) abs(responsets_100.as.slowest)];
% rponsets_eo = [abs(responsets_100.eo.fastest) abs(responsets_100.eo.medfast) abs(responsets_100.eo.slowest)];


%% JAGS code for simple linear regression
model = {
'model {'
    '##########'
    '#Prior of Regression intercept (remaining reaction time)'
    'intercept ~ dnorm(.3, pow(.25,-2))'
    '#Prior of effect of RP latency on reaction time'
    'beta1 ~ dnorm(1,pow(5,-2)) #''Informative'' prior for Bayes Factor calculation'
    '#Prior of Variability in reaction time'
    'tersubsd ~ dgamma(.2,1)'

    '##########'
    '#Linear regression (normal likelihood)'
    'for (i in 1:N) {'
        'rtpercentile[i] ~ dnorm(intercept + RPonset[i]*beta1, pow(tersubsd,-2) )'
    '}'
'}'
};

%% Code for Trinity

% Track these variables
params = {'intercept', 'beta1', 'tersubsd'};

initstruct = @()struct(...
'beta1', rand*4 - 2, ...
'intercept', rand, ...
'tersubsd', rand*.09 + .01);

% Organize data together
dataRTas.rtpercentile = RTmedian_as/1000; %convert to secs
dataRTas.RPonset = rponsets_as/1000; %convert to secs
dataRTas.N = length(dataRTas.rtpercentile);
dataRTeo.rtpercentile = RTmedian_eo/1000; %convert to secs
dataRTeo.RPonset = rponsets_eo/1000; %convert to secs
dataRTeo.N = length(dataRTeo.rtpercentile);
dataDTas.rtpercentile = DTmedian_as/1000; %convert to secs
dataDTas.RPonset = rponsets_as/1000; %convert to secs
dataDTas.N = length(dataDTas.rtpercentile);
dataDTeo.rtpercentile = DTmedian_eo/1000; %convert to secs
dataDTeo.RPonset = rponsets_eo/1000; %convert to secs
dataDTeo.N = length(dataDTeo.rtpercentile);

%% Run RP onset vs reaction time in action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_regress_expanded_RTas';

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
% title('RP onset versus RT Percentiles')
numerator = density_RT(xi==1);
denominator = normpdf(1,1,5);
bf_RTas = numerator/denominator;
intercept_RTas = prctile(chains_RTas.intercept(:),[2.5,50,97.5]);
slope_RTas = prctile(chains_RTas.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the RP vs RT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',intercept_RTas(2),intercept_RTas(1),intercept_RTas(3));
fprintf('The slope parameter of the RP vs RT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',slope_RTas(2),slope_RTas(1),slope_RTas(3));
fprintf('The Bayes Factor of the slope parameter is %.3f ! \n',bf_RTas);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_RTas', 'chains_RTas', 'diagnostics_RTas','info_RTas','params','dataRTas','model','bf_RTas');
system('rm -r wdir');

%% Run RP onset vs decision time in action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_regress_expanded_DTas';

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
% title('RP onset versus RT Percentiles')
numerator = density_DT(xi==1);
denominator = normpdf(1,1,5);
bf_DTas = numerator/denominator;
intercept_DTas = prctile(chains_DTas.intercept(:),[2.5,50,97.5]);
slope_DTas = prctile(chains_DTas.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the RP vs DT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',intercept_DTas(2),intercept_DTas(1),intercept_DTas(3));
fprintf('The slope parameter of the RP vs DT regression in the AS task is %.3f with [%.3f, %.3f] ! \n',slope_DTas(2),slope_DTas(1),slope_DTas(3));
fprintf('The Bayes Factor of slope parameter is %.3f ! \n',bf_DTas);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_DTas', 'chains_DTas', 'diagnostics_DTas','info_DTas','params','dataDTas','model','bf_DTas');
system('rm -r wdir');

%% Run RP onset vs reaction time in execution only (EO)action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_regress_expanded_RTeo';

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
% title('RP onset versus RT Percentiles')
numerator = density_RT(xi==1);
denominator = normpdf(1,1,5);
bf_RTeo = numerator/denominator;
intercept_RTeo = prctile(chains_RTeo.intercept(:),[2.5,50,97.5]);
slope_RTeo = prctile(chains_RTeo.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the RP vs RT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',intercept_RTeo(2),intercept_RTeo(1),intercept_RTeo(3));
fprintf('The slope parameter of the RP vs RT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',slope_RTeo(2),slope_RTeo(1),slope_RTeo(3));
fprintf('The inverse Bayes Factor of the slope parameter is %.3f ! \n',1/bf_RTeo);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_RTeo', 'chains_RTeo', 'diagnostics_RTeo','info_RTeo','params','dataRTeo','model','bf_RTeo');
system('rm -r wdir');

%% Run RP onset vs decision time in execution only (EO) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_regress_expanded_DTeo';

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
% title('RP onset versus RT Percentiles')
numerator = density_DT(xi==1);
denominator = normpdf(1,1,5);
bf_DTeo = numerator/denominator;
intercept_DTeo = prctile(chains_DTeo.intercept(:),[2.5,50,97.5]);
slope_DTeo = prctile(chains_DTeo.beta1(:),[2.5,50,97.5]);
fprintf('The intercept parameter of the RP vs DT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',intercept_DTeo(2),intercept_DTeo(1),intercept_DTeo(3));
fprintf('The slope parameter of the RP vs DT regression in the EO task is %.3f with [%.3f, %.3f] ! \n',slope_DTeo(2),slope_DTeo(1),slope_DTeo(3));
fprintf('The inverse Bayes Factor of the slope parameter is %.3f ! \n',1/bf_DTeo);


save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_DTeo', 'chains_DTeo', 'diagnostics_DTeo','info_DTeo','params','dataDTeo','model','bf_DTeo');
system('rm -r wdir');