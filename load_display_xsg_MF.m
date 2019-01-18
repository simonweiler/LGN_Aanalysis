function load_display_xsg_MF(filterephys,ramp_nr,trace_range)

base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;
redpeak_start       =   100;
redpeak_end         =   349;
bluepeak_start      =   351;
bluepeak_end        =   400;

%filterephys=filter yes or no
%trace_range: range or single value, ex. 1:11; 2:4; or 1, 4 etc. 

directory='F:\LISBOA2019';%change accordingly to  drive 
exp_folder=uipickfiles('FilterSpec',directory);

list=dir([char(exp_folder) '\*.xsg']);%xsg files per cell 
len=length(list);%number of xsg files per cell

            for j=1:len
                load([char(exp_folder) '/' list(j).name],'-mat');%load each xsg file
                iterations(:,j)=header.loopGui.loopGui.iterations;%find out whether mini or ramp recording
            end
            ramp=find(iterations==11);%ramp recordings
            failure1=find(iterations==50);%mini recordings
            failure2=find(iterations==100);%mini recordings
  curr_ramps=ramp(ramp_nr);
      load([char(exp_folder) '/' list(curr_ramps(1)).name],'-mat');%load each xsg file 
      sr = header.ephys.ephys.sampleRate;%check sample rate
      srF = 1/(1000/sr);
      samples_per_sweep = header.ephys.ephys.traceLength*sr;
      timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
      traces=data.ephys.trace_1;%raw ephys trace
      ind_traces=reshape(traces,[length(traces)/11 11]);
 
 for j=1:length(trace_range);
     raw_traces=ind_traces(:,j);
      bs=raw_traces(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
      bs_std=std(bs);%std of baseline trace
      bs_traces=raw_traces-mean(raw_traces(base_start*srF:base_end*srF,:));%subtract baseline
      unfil_traces(:,j)= bs_traces;
      if filterephys==1
      cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
      order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
      type        = 'Butter';
      filt_traces = lowpassfilt(raw_traces, order, cutoff, sr, type);
      bs=filt_traces(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
      bs_std=std(bs);%std of baseline trace
      bs_traces=filt_traces-mean(filt_traces(base_start*srF:base_end*srF,:));%subtract baseline
      final_traces(:,j)=bs_traces;
      end  
   end     

%PLOT
figure;
set(gcf, 'Position', [200, 0, 1000, 500]);
plot(unfil_traces(1:20000,:,:),'Color','k');
hold on;
                y1=get(gca,'ylim');
                x1= redpeak_start*srF;
                hold on;
                p1=plot([x1 x1],y1,'--','Color','r');
                p1.Color(4) = 1;
                hold on;
                y1=get(gca,'ylim');
                x1=redpeak_end*srF;
                hold on;
%                 p2=plot([x1 x1],y1,'--','Color','r');
%                 p2.Color(4) = 1;
%                 hold on;
                %%blue vertical lines
                y1=get(gca,'ylim');
                x1=bluepeak_start*srF;
                hold on;
                p3=plot([x1 x1],y1,'--','Color','b');
                p3.Color(4) = 1;
                hold on;
                y1=get(gca,'ylim');
                x1=bluepeak_end  *srF;
                hold on;
                p4=plot([x1 x1],y1,'--','Color','b');
                p4.Color(4) = 1;
ylabel('Synaptic Input (pA)');
xlabel('Time (ms)');  
figure;
set(gcf, 'Position', [200, 0, 1000, 500]);
plot(final_traces(1:20000,:,:),'Color','k');
hold on;
                y1=get(gca,'ylim');
                x1= redpeak_start*srF;
                hold on;
                p1=plot([x1 x1],y1,'--','Color','r');
                p1.Color(4) = 1;
                hold on;
                y1=get(gca,'ylim');
                x1=redpeak_end*srF;
                hold on;
%                 p2=plot([x1 x1],y1,'--','Color','r');
%                 p2.Color(4) = 0.3;
%                 hold on;
                %%blue vertical lines
                y1=get(gca,'ylim');
                x1=bluepeak_start*srF;
                hold on;
                p3=plot([x1 x1],y1,'--','Color','b');
                p3.Color(4) = 1;
                hold on;
                y1=get(gca,'ylim');
                x1=bluepeak_end  *srF;
                hold on;
                p4=plot([x1 x1],y1,'--','Color','b');
                p4.Color(4) = 1;
ylabel('Synaptic Input (pA)');
xlabel('Time (ms)');
set(gca,'XTick',[0:2000:20000],'XTickLabel',{'0','100','200','300','400','500','600','700','800','900','1000'});
end