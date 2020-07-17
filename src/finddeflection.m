%% this function finds the point of onset of the ERP 

function [onset,erpdiff] = finddeflection(erp,baseline,erpwindow);
maxtime = size(erp,1) - length(baseline);
erpdiff = erp(2:end,:)-erp(1:end-1,:);
for j = 1:size(erp,2)
test = erpdiff(:,j) < 0;
testtime = find(test);
testval = zeros(length(test));
if ~isempty(testtime)
for k = 1:length(testtime)
if testtime(k) < maxtime
testval(testtime(k)) = sum(test(testtime(k)+baseline));
end;
end;
testbaseline = find(testval == length(baseline));
if ~isempty(testbaseline)
onset(j) = erpwindow(testbaseline(1));
else
onset(j) = 0;
end;
else
onset(j) = 0;
end;
end;
onset = onset+length(baseline)/2;