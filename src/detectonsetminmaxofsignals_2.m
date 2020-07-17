%% this script finds the P300 peaks and N200 peak latencies for each participant
%% finds the maximum point of the P300 for each participant in both tasks (Figure 3b)
subjects = {'';'';} %put all the subject IDs in

parchans = [100 110 109 108 128 127 99 119 101 118 139 140 98]; 

for i = 1:length(subjects)
   subid = subjects{i};
   eval(['load ' subid '/' subid '_motorerp_stimlocked.mat']); 
   erp.as.all = (erp.as.left+erp.as.right)/2;
   erp.eo.all = (erp.eo.left+erp.eo.right)/2;
   
   p300maxsignal.as.all = mean(erp.as.all(:,parchans),2);
   p300maxsignal.eo.all = mean(erp.eo.all(:,parchans),2);
   
   [maxsignal.as.all,maxindex.as.all] = max(p300maxsignal.as.all(600:850));
   [maxsignal.eo.all,maxindex.eo.all] = max(p300maxsignal.eo.all(600:850));
   
   p300max.as.all(i) = (maxindex.as.all + 600) - 400;
   p300max.eo.all(i) = (maxindex.eo.all + 600) - 400;
   
end
%% finds the maximum point of P300 when RT is split into tertiles (Figure 6)

for i = 1:length(subjects)
   subid = subjects{i};
   eval(['load ' subid '/' subid '_motorerp_RTthirds_stimlocked.mat']); 
   
   p300maxsignal.as.fastest = mean(erp.as.fastest(:,parchans),2);
   p300maxsignal.eo.fastest = mean(erp.eo.fastest(:,parchans),2);
   p300maxsignal.as.medfast = mean(erp.as.medfast(:,parchans),2);
   p300maxsignal.eo.medfast = mean(erp.eo.medfast(:,parchans),2);
   p300maxsignal.as.slowest = mean(erp.as.slowest(:,parchans),2);
   p300maxsignal.eo.slowest = mean(erp.eo.slowest(:,parchans),2);
   
   [maxsignal.as.fastest,maxindex.as.fastest] = max(p300maxsignal.as.fastest(600:850));
   [maxsignal.eo.fastest,maxindex.eo.fastest] = max(p300maxsignal.eo.fastest(600:850));
   [maxsignal.as.medfast,maxindex.as.medfast] = max(p300maxsignal.as.medfast(600:850));
   [maxsignal.eo.medfast,maxindex.eo.medfast] = max(p300maxsignal.eo.medfast(600:850));
   [maxsignal.as.slowest,maxindex.as.slowest] = max(p300maxsignal.as.slowest(600:850));
   [maxsignal.eo.slowest,maxindex.eo.slowest] = max(p300maxsignal.eo.slowest(600:850));
  
   
   p300max.as.fastest(i) = (maxindex.as.fastest + 600) - 400;
   p300max.eo.fastest(i) = (maxindex.eo.fastest + 600) - 400;
   p300max.as.medfast(i) = (maxindex.as.medfast + 600) - 400;
   p300max.eo.medfast(i) = (maxindex.eo.medfast + 600) - 400;
   p300max.as.slowest(i) = (maxindex.as.slowest + 600) - 400;
   p300max.eo.slowest(i) = (maxindex.eo.slowest + 600) - 400;
   
   
end

%% finds N200 latency between two tasks each subject (Figure 3a)
subjects = {'';''};

n200chans = [107 96 106 97 95 105 160 169 170 159 179 192];

for i = 1:length(subjects)
   subid = subjects{i};
   eval(['load ' subid '/' subid '_erp_N200.mat']); 
   erp.as.all = (erp.as.left+erp.as.right)/2;
   erp.eo.all = (erp.eo.left+erp.eo.right)/2;
   
   n200minsignal.as.all = mean(erp.as.all(:,n200chans),2);
   n200minsignal.eo.all = mean(erp.eo.all(:,n200chans),2);
   
   [minsignal.as.all,minindex.as.all] = min(n200minsignal.as.all(525:675));
   [minsignal.eo.all,minindex.eo.all] = min(n200minsignal.eo.all(525:675));
   
   n200min.as.all(i) = minindex.as.all + 525;
   n200min.eo.all(i) = minindex.eo.all + 525;
   
end

%% finds the N200 latency when RT is split into tertiles

for i = 1:length(subjects)
   subid = subjects{i};
   eval(['load ' subid '/' subid '_N200_RTthirds.mat']); 

   
   n200minsignal.as.fastest = mean(erp.as.fastest(:,n200chans),2);
   n200minsignal.eo.fastest = mean(erp.eo.fastest(:,n200chans),2);
   n200minsignal.as.medfast = mean(erp.as.medfast(:,n200chans),2);
   n200minsignal.eo.medfast = mean(erp.eo.medfast(:,n200chans),2);
   n200minsignal.as.slowest = mean(erp.as.slowest(:,n200chans),2);
   n200minsignal.eo.slowest = mean(erp.eo.slowest(:,n200chans),2);
   
   [minsignal.as.fastest,minindex.as.fastest] = min(n200minsignal.as.fastest(525:675));
   [minsignal.eo.fastest,minindex.eo.fastest] = min(n200minsignal.eo.fastest(525:675));
   [minsignal.as.medfast,minindex.as.medfast] = min(n200minsignal.as.medfast(525:675));
   [minsignal.eo.medfast,minindex.eo.medfast] = min(n200minsignal.eo.medfast(525:675));
   [minsignal.as.slowest,minindex.as.slowest] = min(n200minsignal.as.slowest(525:675));
   [minsignal.eo.slowest,minindex.eo.slowest] = min(n200minsignal.eo.slowest(525:675));
   
   n200min.as.fastest(i) = minindex.as.fastest + 525;
   n200min.eo.fastest(i) = minindex.eo.fastest + 525;
   n200min.as.medfast(i) = minindex.as.medfast + 525;
   n200min.eo.medfast(i) = minindex.eo.medfast + 525;
   n200min.as.slowest(i) = minindex.as.slowest + 525;
   n200min.eo.slowest(i) = minindex.eo.slowest + 525;
end
