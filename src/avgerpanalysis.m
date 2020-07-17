%% finding the onset times of the potentials 

load young_motorerp
baseline = [1:100];
nstd = 2;  %threshhold in terms of standard deviations of baseline.
%collapse first
allerp.eo.all = (allerp.eo.left + allerp.eo.right)/2;
allerp.as.all = (allerp.as.left + allerp.as.right)/2;

%means
avgerp.eo.left = squeeze(mean(allerp.eo.left,1));
avgerp.eo.right = squeeze(mean(allerp.eo.right,1));
avgerp.as.left = squeeze(mean(allerp.as.left,1));
avgerp.as.right = squeeze(mean(allerp.as.right,1));
avgerp.eo.all = squeeze(mean(allerp.eo.all,1));
avgerp.as.all = squeeze(mean(allerp.as.all,1));

%standard deviations
stderp.eo.left = squeeze(std(allerp.eo.left,1));
stderp.eo.right = squeeze(std(allerp.eo.right,1));
stderp.as.left = squeeze(std(allerp.as.left,1));
stderp.as.right = squeeze(std(allerp.as.right,1));
stderp.eo.all = squeeze(std(allerp.eo.all,1));
stderp.as.all = squeeze(std(allerp.as.all,1));
%calculate threshholds 
meanbaseline = mean(avgerp.eo.left(baseline,:),1);
stdbaseline = std(avgerp.eo.left(baseline,:),1);  
cutoffminus.eo.left = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(avgerp.eo.right(baseline,:),1);
stdbaseline = std(avgerp.eo.right(baseline,:),1);  
cutoffminus.eo.right = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(avgerp.as.left(baseline,:),1);
stdbaseline = std(avgerp.as.left(baseline,:),1);  
cutoffminus.as.left = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(avgerp.as.right(baseline,:),1);
stdbaseline = std(avgerp.as.right(baseline,:),1);  
cutoffminus.as.right = -nstd*stdbaseline+meanbaseline;

meanbaseline = mean(avgerp.eo.all(baseline,:),1);
stdbaseline = std(avgerp.eo.all(baseline,:),1);  
cutoffminus.avg.eo.all = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(avgerp.as.all(baseline,:),1);
stdbaseline = std(avgerp.as.all(baseline,:),1);  
cutoffminus.avg.as.all = -nstd*stdbaseline+meanbaseline;


%find deflection points in average
onset.avg.eo.left = finddeflection(avgerp.eo.left,cutoffminus.eo.left,baseline,allerp.window);
onset.avg.eo.right = finddeflection(avgerp.eo.right,cutoffminus.eo.right,baseline,allerp.window);
onset.avg.as.left = finddeflection(avgerp.as.left,cutoffminus.as.left,baseline,allerp.window);
onset.avg.as.right = finddeflection(avgerp.as.right,cutoffminus.as.right,baseline,allerp.window);

onset.avg.as.all = finddeflection(avgerp.as.all,cutoffminus.avg.as.all,baseline,allerp.window);
onset.avg.eo.all = finddeflection(avgerp.eo.all,cutoffminus.avg.eo.all,baseline,allerp.window);

for j = 1:size(allerp.eo.all)
%calculate threshholds for each subject
meanbaseline = mean(squeeze(allerp.eo.left(j,baseline,:)),1);
stdbaseline = std(squeeze(allerp.eo.left(j,baseline,:)),1);  
cutoffminus.eo.left(j,:) = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(squeeze(allerp.eo.right(j,baseline,:)),1);
stdbaseline = std(squeeze(allerp.eo.right(j,baseline,:)),1);  
cutoffminus.eo.right = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(squeeze(allerp.as.left(j,baseline,:)),1);
stdbaseline = std(squeeze(allerp.as.left(j,baseline,:)),1);  
cutoffminus.as.left = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(squeeze(allerp.as.right(j,baseline,:)),1);
stdbaseline = std(squeeze(allerp.as.right(j,baseline,:)),1);  
cutoffminus.as.right = -nstd*stdbaseline+meanbaseline;

meanbaseline = mean(squeeze(allerp.eo.all(j,baseline,:)),1);
stdbaseline = std(squeeze(allerp.eo.all(j,baseline,:)),1);  
cutoffminus.eo.all(j,:) = -nstd*stdbaseline+meanbaseline;
meanbaseline = mean(squeeze(allerp.as.all(j,baseline,:)),1);
stdbaseline = std(squeeze(allerp.as.all(j,baseline,:)),1);  
cutoffminus.as.all(j,:) = -nstd*stdbaseline+meanbaseline;


%find deflection points in average
onset.eo.left(j,:) = finddeflection(squeeze(allerp.eo.left(j,:,:)),cutoffminus.eo.left,baseline,allerp.window);
onset.eo.right(j,:) = finddeflection(squeeze(allerp.eo.right(j,:,:)),cutoffminus.eo.right,baseline,allerp.window);
onset.as.left(j,:) = finddeflection(squeeze(allerp.as.left(j,:,:)),cutoffminus.as.left,baseline,allerp.window);
onset.as.right(j,:) = finddeflection(squeeze(allerp.as.right(j,:,:)),cutoffminus.as.right,baseline,allerp.window);

onset.eo.all(j,:) = finddeflection(squeeze(allerp.eo.all(j,:,:)),cutoffminus.eo.all(j,:),baseline,allerp.window);
onset.as.all(j,:) = finddeflection(squeeze(allerp.as.all(j,:,:)),cutoffminus.as.all(j,:),baseline,allerp.window);
end;
save young_avgmotorerp_2std avgerp stderp onset allerp
