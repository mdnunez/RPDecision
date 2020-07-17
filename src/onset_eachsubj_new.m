%% calculate onset time for resplocked & stim locked for each subject for the RP and then averages across each person's onset
ourchans = [8 9 16 17 23 24 15 186];
subjects = {'';'';};
%% RESP LOCKED -RT thirds (Figure 5a)


for i = 1:length(subjects);
   subid = subjects{i};
   eval(['load ' subid '/' subid '_motorerp_RTthirds_new_resplocked.mat']); 
   %as fastest 
 
   [onset,erpdiff] = finddeflection(mean(erp.as.fastest(201:end,ourchans),2),[1:125],erpwindowfast(201:end));
   responsets.as.fastest(i) = onset; %nanmean(onset(ourchans));

   %as medfast 
   [onset,erpdiff] = finddeflection(mean(erp.as.medfast(201:end,ourchans),2),[1:125],erpwindowmed(201:end));
   responsets.as.medfast(i) = onset;

   %as slowest 
 
   [onset,erpdiff] = finddeflection(mean(erp.as.slowest(201:end,ourchans),2),[1:125],erpwindowslow(201:end));
   responsets.as.slowest(i) = onset;
% 

   %eo fastest 
   [onset,erpdiff] = finddeflection(mean(erp.eo.fastest(201:end,ourchans),2),[1:125],erpwindowfast(201:end));
   responsets.eo.fastest(i) = onset;
% 
   
   %eo medfast 
 
   [onset,erpdiff] = finddeflection(mean(erp.eo.medfast(201:end,ourchans),2),[1:125],erpwindowmed(201:end));
   responsets.eo.medfast(i) = onset;


   %eo slowest
   [onset,erpdiff] = finddeflection(mean(erp.eo.slowest(201:end,ourchans),2),[1:125],erpwindowslow(201:end));
   responsets.eo.slowest(i) = onset;

   %avg onset
   
   respavgonset.as(i) = (responsets.as.fastest(i)+responsets.as.medfast(i) + responsets.as.slowest(i))/3;
   respavgonset.eo(i) = (responsets.eo.fastest(i)+responsets.eo.medfast(i) + responsets.eo.slowest(i))/3;
   
end
%% STIMLOCKED - RT thirds (Figure 5b)
erpwindow = [-400:1500];
for i = 1:length(subjects);
   subid = subjects{i};
   eval(['load ' subid '/' subid '_motorerp_RTthirds_stimlocked.mat']); 
   
   %as fastest
   [onset,erpdiff] = finddeflection(mean(erp.as.fastest(501:end,ourchans),2),[1:100],erpwindow(501:end));
   stimonsets.as.fastest(i) = onset;
   
   %as medfast
   [onset,erpdiff] = finddeflection(mean(erp.as.medfast(501:end,ourchans),2),[1:100],erpwindow(501:end));
   stimonsets.as.medfast(i) = onset;
   %as slowest
   [onset,erpdiff] = finddeflection(mean(erp.as.slowest(501:end,ourchans),2),[1:100],erpwindow(501:end));
   stimonsets.as.slowest(i) = onset;
   
   %eo fastest
   [onset,erpdiff] = finddeflection(mean(erp.eo.fastest(501:end,ourchans),2),[1:100],erpwindow(501:end));
   stimonsets.eo.fastest(i) = onset;
   
   %eo medfast
   [onset,erpdiff] = finddeflection(mean(erp.eo.medfast(501:end,ourchans),2),[1:100],erpwindow(501:end));
   stimonsets.eo.medfast(i) = onset;
   
   %eo slowest
   [onset,erpdiff] = finddeflection(mean(erp.eo.slowest(501:end,ourchans),2),[1:100],erpwindow(501:end));
   stimonsets.eo.slowest(i) = onset;
   
   stimavgonset.as(i) = (stimonsets.as.fastest(i)+stimonsets.as.medfast(i) + stimonsets.as.slowest(i))/3;
   stimavgonset.eo(i) = (stimonsets.eo.fastest(i)+stimonsets.eo.medfast(i) + stimonsets.eo.slowest(i))/3;
end

%% Stimlocked/Resplocked - With Rotation 
%resplocked (Figure 3C)
for i = 1:length(subjects)
   subid = subjects{i};
   eval(['load ' subid '/' subid '_motorerp.mat']); 
   erp.as.all = (erp.as.left+erp.as.right)/2;
   erp.eo.all = (erp.eo.left+erp.eo.right)/2;
    % eo left
    [onset,erpdiff] = finddeflection(mean(erp.eo.left(201:end,ourchans),2),[1:125],erpwindow(201:end));
    responsets.eo.left(i) = onset;
   %eo right
    [onset,erpdiff] = finddeflection(mean(erp.eo.right(201:end,ourchans),2),[1:125],erpwindow(201:end));
    responsets.eo.right(i) = onset;
    %eo all
    [onset,erpdiff] = finddeflection(mean(erp.eo.all(201:end,ourchans),2),[1:125],erpwindow(201:end));
    responsets.eo.all(i) = onset;
    %as left
     [onset,erpdiff] = finddeflection(mean(erp.as.left(201:end,ourchans),2),[1:125],erpwindow(201:end));
    responsets.as.left(i) = onset;
   %as right
    [onset,erpdiff] = finddeflection(mean(erp.as.right(201:end,ourchans),2),[1:125],erpwindow(201:end));
    responsets.as.right(i) = onset;
    %as all
    [onset,erpdiff] = finddeflection(mean(erp.as.all(201:end,ourchans),2),[1:125],erpwindow(201:end));
    responsets.as.all(i) = onset;
end

%stimlocked (Figure 4)
for i = 1:length(subjects)
   subid = subjects{i};
   eval(['load ' subid '/' subid '_motorerp_stimlocked.mat']); 
   erp.as.all = (erp.as.left+erp.as.right)/2;
   erp.eo.all = (erp.eo.left+erp.eo.right)/2;
    %eo left
   [onset,erpdiff] = finddeflection(mean(erp.eo.left(501:end,ourchans),2),[1:100],erp.window(501:end));
   stimonsets.eo.left(i) = onset;
    %eo right
   [onset,erpdiff] = finddeflection(mean(erp.eo.right(501:end,ourchans),2),[1:100],erp.window(501:end));
   stimonsets.eo.right(i) = onset;
%     eo all
   [onset,erpdiff] = finddeflection(mean(erp.eo.all(501:end,ourchans),2),[1:100],erp.window(501:end));
   stimonsets.eo.all(i) = onset;
    %as left
   [onset,erpdiff] = finddeflection(mean(erp.as.left(501:end,ourchans),2),[1:100],erp.window(501:end));
   stimonsets.as.left(i) = onset;
    %as right
   [onset,erpdiff] = finddeflection(mean(erp.as.right(501:end,ourchans),2),[1:100],erp.window(501:end));
   stimonsets.as.right(i) = onset;
%    %as all 
   [onset,erpdiff] = finddeflection(mean(erp.as.all(501:end,ourchans),2),[1:100],erp.window(501:end));
   stimonsets.as.all(i) = onset;
end

save allonset_eachsubj.mat stimavgonset stimonsets  respavgonset responsets 
