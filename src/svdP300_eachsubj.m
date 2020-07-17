%% this script does the svd on the erp that was for the P300 and plots out the components
close all
clear all
subid = '';
saveplotxs = 0;  %set to 1 if you want it to save out your plotxs to a file 
% chanlists = getchan('egi256');  %you can change this to 'cramereeg' to match jessica or 'egi256' to match Bill
eval(['load ' subid '/' subid '_motorerp_RTthirds_stimlocked.mat erp']);
time = [200:450];
colorlim = [-0.1 0.1];
samps = time - erp.window(1)+1;
[uf,sf,vf] = svd(erp.as.fastest(samps,:)); 
[um,sm,vm] = svd(erp.as.medfast(samps,:)); 
[us,ss,vs] = svd(erp.as.slowest(samps,:)); 
HM = load('~/data10/finalheadmodel/EGI256_MNIHEADMODEL');
datatype = 'eeg';
cmap = hotncold;
figure
% us(:,1) = -us(:,1);
% vs(:,1) = -vs(:,1);
% uf(:,1) = -uf(:,1);
% vf(:,1) = -vf(:,1);
% um(:,1) = -um(:,1);
% vm(:,1) = -vm(:,1);

figure
plotx(time,[sqrt(sf(1,1))*uf(:,1) sqrt(sm(1,1))*um(:,1) sqrt(ss(1,1))*us(:,1)]);
title('AS Task')
h = maketopo(vf(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('AS Task: Fastest')
h = maketopo(vm(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('AS Task: Medium');
h = maketopo(vs(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('AS Task: Slowest')
as.uf = sqrt(sf(1,1))*uf;
as.um = sqrt(sm(1,1))*um;
as.us = sqrt(ss(1,1))*us;
[uf,sf,vf] = svd(erp.eo.fastest(samps,:)); 
[um,sm,vm] = svd(erp.eo.medfast(samps,:)); 
[us,ss,vs] = svd(erp.eo.slowest(samps,:)); 
% us(:,1) = -us(:,1);
% vs(:,1) = -vs(:,1);
% uf(:,1) = -uf(:,1);
% vf(:,1) = -vf(:,1);
% um(:,1) = -um(:,1);
% vm(:,1) = -vm(:,1);
figure
plotx(time,[sqrt(sf(1,1))*uf(:,1) sqrt(sm(1,1))*um(:,1) sqrt(ss(1,1))*us(:,1)]);
title('EO Task');
figure

h = maketopo(vf(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('EO Task: Fastest')
h = maketopo(vm(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('EO Task: Medium');
h = maketopo(vs(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('EO Task: Slowest')
eo.uf = sqrt(sf(1,1))*uf;
eo.um = sqrt(sm(1,1))*um;
eo.us = sqrt(ss(1,1))*us;
%%
eval(['load ' subid '/' subid '_motorerp_stimlocked.mat erp']);
erp.eo.all = (erp.eo.left+erp.eo.right)/2;
erp.as.all = (erp.as.left+erp.as.right)/2;
[ua,sa,va] = svd(erp.as.all(samps,:)); 
% ua(:,1) = -ua(:,1);
% va(:,1) = -va(:,1);
figure
plotx(time,sqrt(sa(1,1))*ua(:,1));
title('AS Task');
h = maketopo(va(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('AS Task: All')
as.ua = sqrt(sa(1,1))*ua;


[ua,sa,va] = svd(erp.eo.all(samps,:)); 
% ua(:,1) = -ua(:,1);
% va(:,1) = -va(:,1);
figure
plotx(time,sqrt(sa(1,1))*ua(:,1));
title('EO Task');
h = maketopo(va(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('EO Task: All')
eo.ua = sqrt(sa(1,1))*ua;
