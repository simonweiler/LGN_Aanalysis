function  [neg_failure, pos_failure]=minianalysis(list, idx, pathName, fc, show);
%SW181229
%Function that extracts minis by using the std threshold criterion

%list=      information of folders for each cell SW000XX
%idx=       vector with idx of which recording is mini recording
%pathName=  folder name of cell
%fc=        factor of how many stds the signal should be included 
%show=      show plots or not (1 or 0)

base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;

for i=1:length(idx);
load([char(pathName) '/' list(idx(i)).name],'-mat');
sr = header.ephys.ephys.sampleRate;%check sample rate
srF = 1/(1000/sr);
ephystraces=data.ephys.trace_1;
traces=reshape(ephystraces, 10000, length(ephystraces)/10000);
bs=traces(base_start*srF:base_end*srF,:);
bs_std=std(bs);
bs_traces=bsxfun(@minus, traces, mean(bs));
neg_peak=min(bs_traces(pulse_start*srF:pulse_end*srF+40*srF,:));
pos_peak=max(bs_traces(pulse_start*srF:pulse_end*srF+40*srF,:));
neg_fail=neg_peak<fc*bs_std*(-1);
pos_fail=pos_peak>fc*bs_std;
neg_idx=find(neg_fail==1);
pos_idx=find(pos_fail==1);
neg_m=zeros(1,size(traces,2));
neg_m(neg_idx)=neg_peak(neg_idx);
pos_m=zeros(1,size(traces,2));
pos_m(pos_idx)=pos_peak(pos_idx);

try 
neg_fail(neg_idx)=neg_peak(neg_idx);
pos_fail(pos_idx)=pos_peak(pos_idx);
catch
neg_fail=zeros(length(pos_peak));
pos_fail=zeros(length(pos_peak));
end
neg_failure(:,i)=neg_m;
pos_failure(:,i)=pos_m;

%PLOT
if show==1
figure;
for k=1:size(traces,2)
subplot(5,10,k);
plot(bs_traces(1:2000,k));
end
figure;
plot(neg_peak,'*');
hold on
plot(pos_peak,'*');
hold on; 
plot(repmat(fc*(mean(bs_std))*(-1),length(pos_peak)));
hold on; 
plot(repmat(fc*(mean(bs_std)),length(pos_peak)));
end
end
end
