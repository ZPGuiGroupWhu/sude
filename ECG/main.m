clear;
% Load pre-trained network and ECG data
load ECG_SegNet.mat
data = textread('14046.txt');

% Compute FSST feature
Feature = ECGFeature(data(:,1),250);

% Obtain ECG waves
SegPred = classify(net,Feature);

% Filter out noisy waves
for i = 2:numel(SegPred)-1
    if(SegPred(i-1)==SegPred(i+1)&&SegPred(i+1)~=SegPred(i))
        SegPred(i) = SegPred(i+1);
    end
end

% Obtain ECG segmentation markers
locs = [];
for i = 2:numel(SegPred)
    if(SegPred(i)=='QRS'&&SegPred(i-1)~='QRS')
        locs = [locs;i];
    end
end

% Eliminate incorrect P waves
% 'thres' denotes the threshold of peak amplitude of P wave
thres = 0.4;
rm_id = [];
for i=1:numel(locs)-1
    temp_id = find(SegPred(locs(i):locs(i+1)-1)=='QRS');
    if(length(temp_id) < 3)
        rm_id = [rm_id;i];
    else
        [p_val,~] = findpeaks(data(locs(i)+temp_id-1,1));
        if(isempty(p_val)||max(p_val) < thres)
            rm_id = [rm_id;i];
        end
    end
end
if(length(locs(end):length(data)) < 3)
        rm_id = [rm_id;length(locs)];
else
    [p_val,~] = findpeaks(data(locs(end):end,1));
    if(isempty(p_val)||max(p_val) < thres)
        rm_id = [rm_id;length(locs)];
    end
end
locs(rm_id) = [];

% Plot ECG, segmentation markers, and anomalies
M = signalMask(SegPred);
plotsigroi(M,data(1:3000,1))
hold on;
plot(locs,data(locs,1),'bo');
hold on;
id_n = find(data(:,2)==1);
plot(id_n,data(id_n,1),'go');

% Compute eight ECG features
SegFeature = zeros(length(locs)-1,9);
Key_Pts = zeros(length(locs),1);
for i=1:length(locs)-1
   sub_sig_id = locs(i):locs(i+1)-1;
   sub_sig_val = data(locs(i):locs(i+1)-1,1);
   sub_sig_lab = SegPred(locs(i):locs(i+1)-1);
   SegFeature(i,1) = min(sub_sig_val);
   SegFeature(i,2) = max(sub_sig_val);
   SegFeature(i,7) = length(sub_sig_val);
   [~, peak_loc] = findpeaks(sub_sig_val);
   [~, QRS_id] = max(sub_sig_val);
   T_id = intersect(find(sub_sig_lab=='T'),peak_loc);
   P_id = intersect(find(sub_sig_lab=='P'),peak_loc);
   if (~isempty(T_id)&&~isempty(P_id))
       [~, max_T] = max(sub_sig_val(T_id));
       [~, max_P] = max(sub_sig_val(P_id));
       SegFeature(i,3) = sub_sig_val(T_id(max_T));
       SegFeature(i,4) = sub_sig_val(P_id(max_P));
       SegFeature(i,5) = T_id(max_T)-QRS_id+1;
       SegFeature(i,6) = P_id(max_P)-QRS_id+1;
   end

   Key_Pts(i) = sub_sig_id(QRS_id);
   if(length(find(data(locs(i):locs(i+1)-1,2)))/length(sub_sig_val) > 0.5)
        SegFeature(i,9) = 0;
   else
        SegFeature(i,9) = 1;
   end
end
[~, max_QRS_end] = max(data(locs(end)+1:end,1));
Key_Pts(end) = locs(end) + max_QRS_end;
for i=1:length(locs)-1
    SegFeature(i,8) = Key_Pts(i+1)-Key_Pts(i)+1;
end
% Save features as csv file
csvwrite('ECG_fea.csv',SegFeature);