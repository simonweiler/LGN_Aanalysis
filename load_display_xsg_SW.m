function load_display_xsg_SW(filterephys,trace_range)

%load and display desired trace and filter them: IMPORTANT!!! WORKS ONLY FOR SW data 
%filterephys=filter yes or no
%trace_range: range or single value, ex. 1:11; 2:4; or 1, 4 etc. 

directory='F:\dLGN\example data\';%change accordingly to  drive 
exp_folder=uipickfiles('FilterSpec',directory);
base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;
redpeak_start       =   100;
redpeak_end         =   349;
bluepeak_start      =   351;
bluepeak_end        =   400;

list=dir([char(exp_folder) '\*.xsg']);%xsg files per cell 
len=length(list);%number of xsg files per cell
    
 for j=1:length(trace_range);
      load([char(exp_folder) '/' list(trace_range(j)).name],'-mat');%load each xsg file 
      sr = header.ephys.ephys.sampleRate;%check sample rate
      srF = 1/(1000/sr);
      raw_traces=data.ephys.trace_1;
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
plot(unfil_traces(1:10000,:,:),'Color','k');
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
plot(final_traces(1:10000,:,:),'Color','k');
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
set(gca,'XTick',[0:1000:10000],'XTickLabel',{'0','100','200','300','400','500','600','700','800','900','1000'});
end