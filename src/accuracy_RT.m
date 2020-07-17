%% this saves out the accuracy and response time (RT) for each subject for both tasks 
% figure 2
clear
subid = ''; %subject ID goes here for one subject
eval(['load ' subid '/' subid '_cleaned.mat responseside correctresponses responsetime blockorder']);

%AS = find(segdata.blockorder == 3);
%EO = ~find(segdata.blockorder == 3)

eo = [41:80 121:160 201:240 281:320];
as = setdiff(1:320,eo);
%%Calcuate Accuracy for Action Selection & Execute Only
segdata.expinfo.accuracy = (responseside == correctresponses);
segdata.expinfo.correct = (responseside == correctresponses);
segdata.expinfo.correct;
segdata.expinfo.as_acc = nanmean(segdata.expinfo.correct(as));
segdata.expinfo.eo_acc = nanmean(segdata.expinfo.correct(eo));

%%Calculate RT 
segdata.expinfo.eo_avgrt = nanmean(responsetime(eo));
segdata.expinfo.eo_medrt = nanmedian(responsetime(eo));
segdata.expinfo.eo_rt_right = nanmean(responsetime(find(blockorder == 2)));
segdata.expinfo.eo_rt_left = nanmean(responsetime(find(blockorder == 1)));

segdata.expinfo.as_avgrt = nanmean(responsetime(as));
segdata.expinfo.as_medrt = nanmedian(responsetime(as));
segdata.expinfo.as_rt_right = nanmean(responsetime(intersect(as,find(responseside==2))));
segdata.expinfo.as_rt_left = nanmean(responsetime(intersect(as,find(responseside==1))));

eval(['save -v7.3 ' subid '/' subid '_avgrt.mat segdata' ]);