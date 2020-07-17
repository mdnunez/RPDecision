%puts together median Response Time & onset time(resp_locked) - Figure 8
%% this puts together a file with all the RT and DT median (tertile split) for each subject so you can run the regression model
% this is the same script as regression_rtthirds_p300
% once you have the median_RT_DT_thirds.mat file this can be commented out
subjects = {'ASTK';'CAME';'CHEM';'HEIA';'KHAF';'MAHM';'MATT';'MEJJ';'ORBC';'ORTD';'OTAH';'TRAT';'VEGL';'YAND';};

for i = 1:length(subjects)
    subid = subjects{i};
    eval(['load ' subid '/' subid '_cleaned.mat responsetime']); 
    eval(['load ' subid '/' subid '_RTthirds.mat']); 
    eval(['load ' subid '/' subid '_dt_thirds_NEW.mat']);
    
    RTmedian.as.fastest(i) = nanmedian(responsetime(trials.as.fastest));
    RTmedian.as.medfast(i) = nanmedian(responsetime(trials.as.medfast));
    RTmedian.as.slowest(i) = nanmedian(responsetime(trials.as.slowest));
    RTmedian.eo.fastest(i) = nanmedian(responsetime(trials.eo.fastest));
    RTmedian.eo.medfast(i) = nanmedian(responsetime(trials.eo.medfast));
    RTmedian.eo.slowest(i) = nanmedian(responsetime(trials.eo.slowest));
    
    DTmedian.as.fastest(i) = nanmedian(dt.as.fastest);
    DTmedian.as.medfast(i) = nanmedian(dt.as.medfast);
    DTmedian.as.slowest(i) = nanmedian(dt.as.slowest);
    DTmedian.eo.fastest(i) = nanmedian(dt.eo.fastest);
    DTmedian.eo.medfast(i) = nanmedian(dt.eo.medfast);
    DTmedian.eo.slowest(i) = nanmedian(dt.eo.slowest);



end

save median_RT_DT_thirds_NEWmodel.mat RTmedian DTmedian
%% this script runs the regression model for RT/DT time against RP duration on tertile split; 
% this also spits out the scatter plots, with the R-squared and equation in
% the plot
load('rp_onsets.mat');
load('median_RT_DT_thirds_NEWmodel.mat');
% load('p300maximums.mat');

DT.as.median.all = [DTmedian.as.fastest DTmedian.as.medfast DTmedian.as.slowest];
DT.eo.median.all = [DTmedian.eo.fastest DTmedian.eo.medfast DTmedian.eo.slowest];
RT.as.median.all = [RTmedian.as.fastest RTmedian.as.medfast RTmedian.as.slowest];
RT.eo.median.all = [RTmedian.eo.fastest RTmedian.eo.medfast RTmedian.eo.slowest];
onsets.as.all = [abs(rponset.as.fastest) abs(rponset.as.medfast) abs(rponset.as.slowest)];
onsets.eo.all = [abs(rponset.eo.fastest) abs(rponset.eo.medfast) abs(rponset.eo.slowest)];


LM.as.DT = fitlm(onsets.as.all,DT.as.median.all);
LM.eo.DT = fitlm(onsets.eo.all,DT.eo.median.all);

LM.as.RT = fitlm(onsets.as.all,RT.as.median.all);
LM.eo.RT = fitlm(onsets.eo.all,RT.eo.median.all);

beta.as.DT = LM.as.DT.Coefficients.Estimate;
beta.eo.DT = LM.eo.DT.Coefficients.Estimate;
beta.as.RT = LM.as.RT.Coefficients.Estimate;
beta.eo.RT = LM.eo.RT.Coefficients.Estimate;

range.as = [min(onsets.as.all) max(onsets.as.all)];
range.eo = [min(onsets.eo.all) max(onsets.eo.all)];

model.as.DT = beta.as.DT(1)+range.as*beta.as.DT(2);
model.eo.DT = beta.eo.DT(1)+range.eo*beta.eo.DT(2);
model.as.RT = beta.as.RT(1)+range.as*beta.as.RT(2);
model.eo.RT = beta.eo.RT(1)+range.eo*beta.eo.RT(2);

figure
h = plot(abs(rponset.as.fastest),DTmedian.as.fastest,'ko');
set(h,'MarkerFaceColor','r')
set(h,'MarkerSize',10)
set(gca,'FontSize',12)
hold on
h = plot(abs(rponset.as.medfast),DTmedian.as.medfast,'ko');
set(h,'MarkerFaceColor','b')
set(h,'Marker','d');
set(h,'MarkerSize',10)
hold on
h = plot(abs(rponset.as.slowest),DTmedian.as.slowest,'ko');
set(h,'MarkerFaceColor','g')
set(h,'Marker','s');
set(h,'MarkerSize',10)
ylabel('Decision Time')
xlabel('RP Duration')
plot([0 1600],[0 1600],'k-.')
h = plot(range.as,model.as.DT,'k-');
set(h,'LineWidth',1)
axis([0 1600 0 1600])
legend('AS:Fastest','AS:Middle','AS:Slowest','Location','Northwest')
hold on
thestring = 'DT = %.f + %.2f * RP Duration';
thestring = sprintf(thestring,beta.as.DT(1),beta.as.DT(2));
text(50,1200,thestring,'Fontsize',12);
thestring = 'R^2 = %.2f';
thestring = sprintf(thestring, LM.as.DT.Rsquared.Ordinary);
text(50,1125,thestring,'Fontsize',12);



figure
h = plot(abs(rponset.eo.fastest),DTmedian.eo.fastest,'ko');
set(h,'MarkerFaceColor','r')
set(h,'MarkerSize',10)
set(gca,'FontSize',12)
hold on
h = plot(abs(rponset.eo.medfast),DTmedian.eo.medfast,'ko');
set(h,'MarkerFaceColor','b')
set(h,'Marker','d');
set(h,'MarkerSize',10)
hold on
h = plot(abs(rponset.eo.slowest),DTmedian.eo.slowest,'ko');
set(h,'MarkerFaceColor','g')
set(h,'Marker','s');
set(h,'MarkerSize',10)
ylabel('Decision Time')
xlabel('RP Duration')
plot([0 1600],[0 1600],'k-.')
h = plot(range.eo,model.eo.DT,'k-');
set(h,'LineWidth',1)
axis([0 1600 0 1600])
legend('EO:Fastest','EO:Middle','EO:Slowest','Location','NorthWest')
thestring = 'DT = %.f + %.2f * RP Duration';
thestring = sprintf(thestring,beta.eo.DT(1),beta.eo.DT(2));
text(50,1200,thestring,'Fontsize',12);
thestring = 'R^2 = %.2f';
thestring = sprintf(thestring, LM.eo.DT.Rsquared.Ordinary);
text(50,1125,thestring,'Fontsize',12);

figure
h = plot(abs(rponset.as.fastest),RTmedian.as.fastest,'ko');
set(h,'MarkerFaceColor','r')
set(h,'MarkerSize',10)
set(gca,'FontSize',12)
hold on
h = plot(abs(rponset.as.medfast),RTmedian.as.medfast,'ko');
set(h,'MarkerFaceColor','b')
set(h,'MarkerSize',10)
set(h,'Marker','d');
hold on
h = plot(abs(rponset.as.slowest),RTmedian.as.slowest,'ko');
set(h,'MarkerFaceColor','g')
set(h,'MarkerSize',10)
set(h,'Marker','s');
ylabel('Response Time')
xlabel('RP Duration')
plot([0 1600],[0 1600],'k-.')
h = plot(range.as,model.as.RT,'k-');
set(h,'LineWidth',1)
axis([0 1600 0 1600])
legend('AS:Fastest','AS:Middle','AS:Slowest','Location','NorthWest')
thestring = 'RT = %.f + %.2f * RP Duration';
thestring = sprintf(thestring,beta.as.RT(1),beta.as.RT(2));
text(50,1200,thestring,'Fontsize',12);
thestring = 'R^2 = %.2f';
thestring = sprintf(thestring, LM.as.RT.Rsquared.Ordinary);
text(50,1125,thestring,'Fontsize',12);

figure
h = plot(abs(rponset.eo.fastest),RTmedian.eo.fastest,'ko');
set(h,'MarkerFaceColor','r')
set(h,'MarkerSize',10)
set(gca,'FontSize',12)
hold on
h = plot(abs(rponset.eo.medfast),RTmedian.eo.medfast,'ko');
set(h,'MarkerFaceColor','b')
set(h,'MarkerSize',10)
set(h,'Marker','d');
hold on
h = plot(abs(rponset.eo.slowest),RTmedian.eo.slowest,'ko');
set(h,'MarkerFaceColor','g')
set(h,'MarkerSize',10)
set(h,'Marker','s');
ylabel('Response Time')
xlabel('RP Duration')
plot([0 1600],[0 1600],'k-.')
h = plot(range.eo,model.eo.RT,'k-');
set(h,'LineWidth',1)
axis([0 1600 0 1600])
legend('EO:Fastest','EO:Middle','EO:Slowest','Location','NorthWest')
thestring = 'RT = %.f + %.2f * RP Duration';
thestring = sprintf(thestring,beta.eo.RT(1),beta.eo.RT(2));
text(50,1200,thestring,'Fontsize',12);
thestring = 'R^2 = %.2f';
thestring = sprintf(thestring, LM.eo.RT.Rsquared.Ordinary);
text(50,1125,thestring,'Fontsize',12);
