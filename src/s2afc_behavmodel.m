%S2AFC_BEHAVMODEL - Finds posterior distributions of parameters in a Hierarchical Drift Diffusion Model in JAGS with a lapse process

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  5/23/17        Michael Nunez           Converted from s2afc_model1.m
%  8/17/18        Michael Nunez               Removal of older participants
%                                        Addition of lapse process


%% Initial
if exist('~/data10') == 7
    data10loc = '~/data10';
elseif exist('/data10') == 7
    data10loc = '/data10';
end
youngloc = [data10loc,'/Stroke2AFC/young'];
jagsinloc = [data10loc,'/Stroke2AFC/summarydata/youngbehavdata.mat'];

saveloc = [data10loc,'/Stroke2AFC/jagsout/'];

% %% Prepare input
% if ~(exist(jagsinloc)==2)
%     fprintf('Extracting data for JAGS...\n');
%     acc = []; %accuracies
%     rt = []; %reaction times (secs)
%     sub = []; %subject index
%     subname = {}; %subject designation
%     aseo = []; %1 - AS, 2 - EO
%     responseside = []; %side of response
%     train =[]; %1 -training data, 0 - test data

%     subtrack = 0;
%         dirout = dir(youngloc);
%         for n=1:numel(dirout)
%             if (dirout(n).isdir == 1) & (sum(ismember(dirout(n).name,'A':'Z')) == 4)
%                 if exist([youngloc,sprintf('/%s/%s_cleaned.mat',dirout(n).name,dirout(n).name)]) == 2
%                     subtrack = subtrack + 1;
%                     subname{subtrack} = dirout(n).name;
%                     fprintf('Loading data for subject %s...\n',subname{subtrack});
%                     eeg = load([youngloc,sprintf('/%s/%s_cleaned.mat',...
%                         subname{subtrack},subname{subtrack})]);
%                     responseside = [responseside, eeg.responseside];
%                     acc = [acc (eeg.responseside == eeg.correctresponses)];
%                     rt = [rt eeg.responsetime];
%                     datalen = length(eeg.responsetime);
%                     sub = [sub subtrack*ones(1,datalen)];
%                     tempaseo = ones(1,datalen);
%                     tempaseo([41:80 121:160 201:240 281:datalen]) = 2; %Note that this may need to be changed in the future
%                     aseo = [aseo tempaseo];
%                     rng(sum(subname{subtrack})); %Set random seed that it consistent for each subject (uses subject names)
%                     testindx = randperm(datalen);
%                     temptrain = ones(1,datalen);
%                     temptrain(testindx(1:round(datalen*.33))) = 0;
%                     train = [train temptrain];
%                 end
%             end
%         end
%     fprintf('Saving extracted data...\n');
%     save(jagsinloc,'acc','rt','sub','subname','aseo','train');
% end

load(jagsinloc);

train = logical(train);
newtrain = train & ~isnan(rt); %Remove NANs
data.y = (2*acc - 1).*(rt/1000);
data.subject = sub;
data.condition = aseo;
data.nconds = 2;
data.nsubs = max(data.subject);

data.minrt = zeros(2,data.nsubs);
data.maxrt = zeros(2,data.nsubs);
for c=1:2,
    for n=1:data.nsubs,
        data.minrt(c,n) = nanmin(abs(data.y(aseo==c & data.subject == n)));
        data.maxrt(c,n) = nanmax(abs(data.y(aseo==c & data.subject == n)));
    end
end

data.condition = data.condition(newtrain);
data.subject = data.subject(newtrain);
data.y = data.y(newtrain);
data.N = length(data.y);

%% JAGS code for the diffusion model
% 3 Parameter Model with the effect of motorERP latency on all parameters across subjects, split by condition
model = {
'model {'
    '##########'
    '#Fixed Parameters'
    'beta <- .5'
    '##########'
    '#Between-subject variability in non-decision time'
    'tersubsd ~ dgamma(.2,1)'

    '#Between-subject variability in drift'
    'deltasubsd ~ dgamma(1,1)'

    '#Between-subject variability in boundary separation'
    'alphasubsd ~ dgamma(1,1)'

    '##########'
    '#Block-level parameters'
    '##########'
    'for (k in 1:nconds) {'

        '#Condition-level non-decision time'
        'tercond[k] ~ dnorm(.3, pow(.25,-2))'

        '#Condition-level drift rate'
        'deltacond[k] ~ dnorm(1, pow(2, -2))'

        '#Condition-level boundary separation'
        'alphacond[k] ~ dnorm(1, pow(.5,-2))'

        '#Subject-level parameters'
        'for (sub in 1:nsubs) {'
            '#Subject-level non-decision time'
            'tersub[k,sub] ~ dnorm(tercond[k],'
               'pow(tersubsd, -2))T(0,1)'

            '#Subject-level drift rate'
            'deltasub[k,sub] ~ dnorm(deltacond[k],'
                'pow(deltasubsd, -2))T(-9, 9)'

            '#Subject-level boundary separation'
            'alphasub[k,sub] ~ dnorm(alphacond[k],'
                'pow(alphasubsd, -2))T(.1,5)'

            '#EEGsession-level lapse trials'
            'probsub[k, sub, 1:2] ~ ddirch(c(1,1))'
        '}'
    '}'
    '##########'
    '#Wiener likelihoods'
    'for (i in 1:N) {'
        '# Log density for DDM process'
        'ldcomp[i, 1] <- dlogwiener(y[i],alphasub[condition[i],subject[i]],'
        'tersub[condition[i],subject[i]],'
        'beta,'
        'deltasub[condition[i],subject[i]])'

        '# Log density for lapse trials (negative max RT to positive max RT)'
        'ldcomp[i, 2] <- logdensity.unif(y[i], -maxrt[condition[i],subject[i]], maxrt[condition[i],subject[i]])'

        '# Select one of these two densities (Mixture of nonlapse and lapse trials)'
        'density[i] <- exp(ldcomp[i, componentchosen[i]] - Constant)'
        
        '# Generate a likelihood for the MCMC sampler using a trick to maximize density value'
        'Ones[i] ~ dbern(density[i])'

        '# Probability of mind wandering trials (lapse trials)'
        'componentchosen[i] ~ dcat(probsub[condition[i],subject[i],1:2])'
    '}'
'}'
};

%% Code for Trinity

data.Constant = 10;
data.Ones = ones(1,data.N);

% Track these variables
params = {'alphasub', 'deltasub', 'tersub', 'probsub',...
             'alphacond', 'deltacond', 'tercond', ...
             'alphasubsd', 'deltasubsd', 'tersubsd'};

initstruct = @()struct(...
'alphasub', rand(data.nconds,data.nsubs)*.5 + 1.5, ...
'deltasub', rand(data.nconds,data.nsubs)*8 -4, ...
'tersub', rand(data.nconds,data.nsubs)*.3, ...
'alphacond', rand(1,data.nconds)*.5 + 1.5, ...
'deltacond', rand(1,data.nconds)*8 -4, ...
'tercond', rand(1,data.nconds)*.3, ...
'alphasubsd', rand*.9 + .1,...
'deltasubsd', rand*2.9 + .1,...
'tersubsd', rand*.09 + .01,...
'probsub', cat(3,zeros(data.nconds,data.nsubs),ones(data.nconds,data.nsubs)));

%% Create model name

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modeltype = 'basic_lapse';

modelname = [modeltype,timestr];

%% Run JAGS

fprintf('Building JAGS model %s and saving output...\n',modelname);

nsamples = 5e3;
nburnin = 2e3;
nchains =3;
thin =10;
verbosity =1;
parallelit = 1; %Set this to 1 if GNU Parallel is installed
maxcores = 3;
modules = {'wiener' 'dic'};

tic
[stats, chains, diagnostics, info] = callbayes('jags', ...
    'model', model, ...
    'data', data, ...
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

save(sprintf('%sjagsmodel%s.mat',saveloc,modelname),'stats', 'chains', 'diagnostics','info','params','data');

