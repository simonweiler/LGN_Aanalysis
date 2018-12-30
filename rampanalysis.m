function  [blue_ramp, red_ramp]=rampanalysis(list, idx, pathName, fc, show, ramp_rtrace);
%SW181229
%function to extract synaptic current peak, integral and photodiode signal
%for blue and red laser

%fct inputs
%   list= list of xsg files per cell 
%   idx=  indices for ramp recodings
%   pathName=folder name of cell 
%   fc= factor for threshold std 
%   show= display plot or not (1 or 0)
%   ramp_rtrace= extract and save raw traces or not (1 or 0)

%define temporal windows  
base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;

redpeak_start       =   100;
redpeak_end         =   349;
bluepeak_start      =   351; 
bluepeak_end        =   400; 


%create vector with start and end point for each ramp within the cell recording
runramp=1:11:length(idx);
runramp=[runramp runramp(end)+11];
%plot if wanted
if show==1
figure;
set(gcf, 'Position', [200, 0, 1500, 1000]);
end
%load each ramp per cell consecutively and extract relevant values such as
%snaptic current peak, integral and photodiode signal for blue and red
%temporal windows
for j=1:(length(idx)/11)% how many ramps in total; loop across ramps per cell 
    counter=1;
    if show==1
    subplot(2,(length(idx)/11)-2,j);
    end
    for i=runramp(j):runramp(j+1)-1;%within each ramp load xsg files (11 in total per ramp)
    load([char(pathName) '/' list(idx(i)).name],'-mat');
    sr = header.ephys.ephys.sampleRate;%check sample rate
    srF = 1/(1000/sr);
    traces=data.ephys.trace_1;%raw ephys trace
    photodiode=data.acquirer.trace_1;%photodiode (PD) signal
    blue_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{4,counter+1}.amplitude;%blue laser amplitude set in ephus 
    try
    red_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{2,counter+1}.amplitude;%red laser amplitude set in ephus
    catch 
    red_amp(j,counter)=0;
    end
    bs=traces(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
    bs_std=std(bs);%std of baseline trace
    bs_traces=traces-mean(traces(base_start*srF:base_end*srF,:));%subtract baseline 
    
    %for first window
    neg_peak1(j,counter)=min(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%negative peak within the red stimulation window 
    pos_peak1(j,counter)=max(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%positive peak within the red stimulation window 
    integ1(j,counter)=trapz(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%Integral within the red stimulation window 
    neg_fail1(j,counter)=neg_peak1(j,counter)<fc*bs_std*(-1);%vector with binary values when neg peaks crossed definded std threshold
    pos_fail1(j,counter)=pos_peak1(j,counter)>fc*bs_std;%vector with binary values when pos peaks crossed definded std threshold
    %photodiode 
    PD1(j,counter)=max(photodiode(redpeak_start*srF:redpeak_end*srF,:));%max values of PD signal within the red stimulation window  
    
    %for second window (same extraction as above for blue laser window 
    neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
    pos_peak2(j,counter)=max(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
    integ2(j,counter)=trapz(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
    neg_fail2(j,counter)=neg_peak2(j,counter)<fc*bs_std*(-1);
    pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
    %photodiode 
    PD2(j,counter)=max(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
    %ephys_traces
    ephys_traces(j,counter,:)=bs_traces;
    
    counter=counter+1;
    traces=[];
    %%%%%%%%%%%%%%plot
    if show==1
    plot(bs_traces(1:10000,:),'linewidth',1,'Color',[0 0 0]+0.05*counter);
    hold on;
    ylabel('Synaptic input (pA)');
    xlabel('Time (ms)');
    end
    end
    if show==1;
    %%red vertical lines
    hold on;
    y1=get(gca,'ylim');
    x1= redpeak_start*srF;
    hold on;
    p1=plot([x1 x1],y1,'--','Color','r');
    p1.Color(4) = 0.3;
    hold on;
    y1=get(gca,'ylim');
    x1=redpeak_end*srF;
    hold on;
    p2=plot([x1 x1],y1,'--','Color','r');
    p2.Color(4) = 0.3;
    hold on;
    %%blue vertical lines
    y1=get(gca,'ylim');
    x1=bluepeak_start*srF;
    hold on;
    p3=plot([x1 x1],y1,'--','Color','b');
    p3.Color(4) = 0.3;
    hold on;
    y1=get(gca,'ylim');
    x1=bluepeak_end  *srF;
    hold on;
    p4=plot([x1 x1],y1,'--','Color','b');
    p4.Color(4) = 0.3;
    end
end
 %%%%%%%%%%%%%%%%%%%%%  output %%%%%%%%%%%%%%%%%
red_ramp.neg_peak1=neg_peak1;
red_ramp.pos_peak1=pos_peak1;
red_ramp.integ1=integ1;
red_ramp.neg_fail1=neg_fail1;
red_ramp.pos_fail1=pos_fail1;
red_ramp.PD=PD1;
red_ramp.laser_amp=red_amp;
if ramp_rtrace==1;
red_ramp.ephys_traces=ephys_traces;
end

%create structure with extracted parameters 
blue_ramp.neg_peak2=neg_peak2;
blue_ramp.pos_peak2=pos_peak2;
blue_ramp.integ2=integ2;
blue_ramp.neg_fail2=neg_fail2;
blue_ramp.pos_fail2=pos_fail2;
blue_ramp.PD=PD2;
blue_ramp.laser_amp=blue_amp;

end

