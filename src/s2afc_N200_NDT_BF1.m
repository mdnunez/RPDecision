%s2afc_N200_NDT_BF1.m - Calculates Bayes Factors for simple regression models with N200 as IV and NDT as DV
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
%  10/20/20        Michael Nunez             Adapted from s2afc_BFregression_rtthirds.m

%% Initial

if exist('~/data10') == 7
    data10loc = '~/data10';
elseif exist('/data10') == 7
    data10loc = '/data10';
end

saveloc = [data10loc,'/Stroke2AFC/jagsout/'];

table = readtable('../data/NDT_ERP.csv');


%% JAGS code for simple linear regression with participants as observations
model = {
'model {'
    '##########'
    '#Prior of Regression intercept (remaining reaction time)'
    'intercept ~ dnorm(.3, pow(.25,-2))'
    '#Prior of effect of N200 latency on non-decision time'
    'beta1 ~ dnorm(1,pow(3,-2)) #''Informative'' prior for Bayes Factor calculation'
    '#Prior of Variability in reaction time'
    'tersubsd ~ dgamma(.2,1)'

    '##########'
    '#Linear regression (normal likelihood)'
    'for (i in 1:N) {'
        'ndt[i] ~ dnorm(intercept + N200latency[i]*beta1, pow(tersubsd,-2) )'
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
dataAS.ndt = table.AS_ndt;
dataAS.N200latency = table.AS_n200/1000; %convert to secs
dataAS.N = length(dataAS.ndt);
dataEO.ndt = table.EO_ndt;
dataEO.N200latency = table.EO_n200/1000; %convert to secs
dataEO.N = length(dataEO.ndt);

%% Run N200 latency vs NDT in action selection (AS) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_N200_NDT_AS';

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
[stats_AS, chains_AS, diagnostics_AS, info_AS] = callbayes('jags', ...
    'model', model, ...
    'data', dataAS, ...
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
[density_RT] = ksdensity(chains_AS.beta1(:),xi);
% figure;
% plot(xi,density_RT,'r');
% line(xi,normpdf(xi,1,3));
% title('MRP onset versus RT Percentiles')
numerator = density_RT(xi==1);
denominator = normpdf(1,1,3);
bf_AS = numerator/denominator;
fprintf('The Bayes Factor of the N200 latency versus NDT regression is %.3f in the AS task! \n',bf_AS);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_AS', 'chains_AS', 'diagnostics_AS','info_AS','params','dataAS','model','bf_AS');
system('rm -r wdir');

%% Run N200 latency vs non-decision time in execution only (EO) condition model

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'BF1_N200_NDT_EO';

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
[stats_EO, chains_EO, diagnostics_EO, info_EO] = callbayes('jags', ...
    'model', model, ...
    'data', dataEO, ...
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
[density_DT] = ksdensity(chains_EO.beta1(:),xi);
% figure;
% plot(xi,density_DT,'r');
% line(xi,normpdf(xi,1,3));
% title('MRP onset versus RT Percentiles')
numerator = density_DT(xi==1);
denominator = normpdf(1,1,3);
bf_EO = numerator/denominator;
fprintf('The Bayes Factor of the N200 latency versus NDT regression is %.3f in the EO task! \n',bf_EO);

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats_EO', 'chains_EO', 'diagnostics_EO','info_EO','params','dataEO','model','bf_EO');
system('rm -r wdir');

