%% this script does svd on the erp that was for the N200 and plots out the components for each participant
close all
clear all
subid = '';
saveplotxs = 0;  %set to 1 if you want it to save out your plotxs to a file 
% chanlists = getchan('egi256');  %you can change this to 'cramereeg' to match jessica or 'egi256' to match Bill
eval(['load ' subid '/' subid '_erp_N200.mat erp']);
time = [100:250];
colorlim = [-0.1 0.1];
samps = time - erp.window(1)+1;
HM = load('~/data10/finalheadmodel/EGI256_MNIHEADMODEL');
datatype = 'eeg';
cmap = hotncold;

erp.eo.all = (erp.eo.left+erp.eo.right)/2;
erp.as.all = (erp.as.left+erp.as.right)/2;
[ua,sa,va] = svd(erp.as.all(samps,:)); 
% ua(:,3) = -ua(:,3);
% va(:,3) = -va(:,3);
% ua(:,1) = -ua(:,1);
% va(:,1) = -va(:,1);
% ua(:,2) = -ua(:,2);
% va(:,2) = -va(:,2);
plotx(time,sqrt(sa(1,1))*ua(:,1));
title('AS Task');
h = maketopo(va(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('AS Task: All')
as.ua = sqrt(sa(1,1))*ua;


[ua,sa,va] = svd(erp.eo.all(samps,:)); 
 figure
% ua(:,1) = -ua(:,1);
% va(:,1) = -va(:,1);
% ua(:,2) = -ua(:,2);
% va(:,2) = -va(:,2);
% ua(:,3) = -ua(:,3);
% va(:,3) = -va(:,3);
plotx(time,sqrt(sa(1,1))*ua(:,1));
title('EO Task')
h = maketopo(va(:,1),HM.Model,HM.Electrode,datatype,[],cmap,colorlim);
title('EO Task: All')
eo.ua = sqrt(sa(1,1))*ua;