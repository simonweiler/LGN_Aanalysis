function  [neg_failure, pos_failure PD1 PD2]=minianalysis(list, idx, pathName, fc, show, user);
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
pulse2_start        =   351; 
pulse2_end          =   360; 

for i=1:length(idx);
load([char(pathName) '/' list(idx(i)).name],'-mat');
sr = header.ephys.ephys.sampleRate;%check sample rate
srF = 1/(1000/sr);
ephystraces=data.ephys.trace_1;
if user==0%SW
traces=reshape(ephystraces, 10000, length(ephystraces)/10000);
photodiode=data.acquirer.trace_1;
photodiode=reshape(photodiode, 10000, length(ephystraces)/10000);
else %MF
traces=reshape(ephystraces, 20000, length(ephystraces)/20000);
photodiode=data.acquirer.trace_1;
photodiode=reshape(photodiode, 20000, length(ephystraces)/20000);
end
bs=traces(base_start*srF:base_end*srF,:);
bs_std=std(bs);
bs_traces=bsxfun(@minus, traces, mean(bs));
bs_photodiode=bsxfun(@minus, photodiode, mean(photodiode(base_start*srF:base_end*srF,:)));
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
PD1(:,i)=max(bs_photodiode(pulse_start *srF:pulse_end*srF,:));
PD2(:,i)=max(bs_photodiode(pulse2_start *srF:pulse2_end*srF,:));
%PLOT
if show==1
figure;
for k=1:size(traces,2)
subplot(size(bs_traces,2)/10,10,k);
if user==0%SW
plot(bs_traces(1:2000,k));
else%MF
plot(bs_traces(1:4000,k));
end
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
