%% this function zeros out the the samples you want to use as your baseline for ERP
function newerp = baselinecorrect(erpdata,baselinesamps);
%sets a baseline of 0 in the interval specified in samples 
%INPUT: erpdata - time by channels 
%	baselinesames: samples to use as baseline, e.g., [1:100];
%OUTPUT: newerp - baseline mean has been set to zero 

baseline = mean(erpdata(baselinesamps,:),1);
newerp = erpdata - ones(size(erpdata,1),1)*baseline;
