
%% 
concatenate_ramp_display==1
display_ramp=0

if display_ramp==0
tr=19;
iter=4;
pulse=4;

base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;
redpeak_start       =   100;
redpeak_end         =   349;
bluepeak_start      =   351;
bluepeak_end        =   400;

raw_traces=data{tr,5}.ephys_traces(:,iter,pulse);
      cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
      order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
      type        = 'Butter';
      sr=10000;
      srF=10;
      filt_traces = lowpassfilt(raw_traces, order, cutoff, sr, type);
      
figure;
%set(gcf, 'Position', [200, 0, 1000, 500]);
plot(filt_traces(1:10000,:),'Color','k');
%hold on;plot(data{tr,5}.ephys_traces(1:10000,8:11,4)+50,'Color','k');
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
    p2=plot([x1 x1],y1,'--','Color','b');
    p2.Color(4) = 1;
    hold on;
%     y1=get(gca,'ylim');
%     x1=bluepeak_start*srF;
%     hold on;
%     p3=plot([x1 x1],y1,'--','Color','b');
%     p3.Color(4) = 0.3;
%     hold on;
    y1=get(gca,'ylim');
    x1=bluepeak_end  *srF;
    hold on;
    p4=plot([x1 x1],y1,'--','Color','b');
    p4.Color(4) = 1;
    xlabel('Time (ms)');
    ylabel('Synaptic input (pA)');
    set(gca,'YTick',[-2000:500:500]);
    set(gca,'XTick',[1000:2000:10000],'XTickLabel',{'100','300','500','700','900','1000'});
    set(gca,'FontSize',15);
    set(gca, 'box', 'off');
    if pulse==3 | pulse==4
        legend('NMDA');
        legend boxoff;
         text(7000, min(filt_traces(1:10000,:))-10, ['irradiance red = ' num2str(data{tr,5}.irr_red(pulse,iter)) ' ' 'mW/mm^{2}'],'Color','r');hold on;
         text(7000, min(filt_traces(1:10000,:))-60, ['irradiance blue = ' num2str(data{tr,4}.irr_blue(pulse,iter)) ' ' 'mW/mm^{2}'],'Color','b');
    else
    legend('AMPA');
    legend boxoff;   
    text(7000, min(filt_traces(1:10000,:))+100, ['irradiance red = ' num2str(data{tr,5}.irr_red(pulse,iter)) ' ' 'mW/mm^{2}'],'Color','r');hold on;
    text(7000, min(filt_traces(1:10000,:))+40, ['irradiance blue = ' num2str(data{tr,4}.irr_blue(pulse,iter)) ' ' 'mW/mm^{2}'],'Color','b');
    end
%% 

else
    tr=3;

pulse=2;

base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;
redpeak_start       =   100;
redpeak_end         =   349;
bluepeak_start      =   351;
bluepeak_end        =   400;
    
figure;
set(gcf, 'Position', [200, 0, 1000, 500]);
hold on;
for i=1:11
 raw_traces=data{tr,5}.ephys_traces(:,i,pulse);
      cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
      order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
      type        = 'Butter';
      sr=10000;
      srF=10;
      filt_traces = lowpassfilt(raw_traces, order, cutoff, sr, type);
     

plot(filt_traces(1:10000,:),'Color',[1 1 1]-0.08*i);
%hold on;plot(data{tr,5}.ephys_traces(1:10000,8:11,4)+50,'Color','k');
 
end
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
    p2=plot([x1 x1],y1,'--','Color','b');
    p2.Color(4) = 1;
    hold on;
%     y1=get(gca,'ylim');
%     x1=bluepeak_start*srF;
%     hold on;
%     p3=plot([x1 x1],y1,'--','Color','b');
%     p3.Color(4) = 0.3;
%     hold on;
    y1=get(gca,'ylim');
    x1=bluepeak_end  *srF;
    hold on;
    p4=plot([x1 x1],y1,'--','Color','b');
    p4.Color(4) = 1;
    xlabel('Time (ms)');
    ylabel('Synaptic input (pA)');
    set(gca,'XTick',[0:1000:10000],'XTickLabel',{'0','100','200','300','400','500','600','700','800','900','1000'});
end   
%% CONCATENATE RAMPS

if concatenate_ramp_display==1
    
tr=5;
pulse=2;

base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;
redpeak_start       =   100;
redpeak_end         =   349;
bluepeak_start      =   351;
bluepeak_end        =   400;
    
figure;
set(gcf, 'Position', [200, 0, 1000, 500]);
hold on;
for i=1:11
 raw_traces=data{tr,5}.ephys_traces(:,i,pulse);
      cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
      order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
      type        = 'Butter';
      sr=10000;
      srF=10;
      filt_traces(:,i) = lowpassfilt(raw_traces, order, cutoff, sr, type);
      power_red(:,i)=data{tr,5}.irr_red(pulse,i);
      power_blue(:,i)=data{tr,4}.irr_blue(pulse,i);
     
end
temp=filt_traces(1:10000,:);
concatenate_ramp=horzcat(temp(:));
temp=[];

figure;plot(concatenate_ramp,'Color','k');
 xlabel('Ramp Nr.');
    ylabel('Synaptic input (pA)');
    set(gca,'XTick',[0:10000:110000]);
hold on
red_mov=[1000:10000:110000];
blue_mov=[3500:10000:110000];
for i=1:11;
     y1=get(gca,'ylim');
    x1= red_mov(i);
    hold on;
    p1=plot([x1 x1],y1,'--','Color','r');
    p1.Color(4) = 1;
end

for i=1:11;
     y1=get(gca,'ylim');
    x1= blue_mov(i);
    hold on;
    p1=plot([x1 x1],y1,'--','Color','b');
    p1.Color(4) = 1;
end
in_val=interp1(1:6,power_red(1:6),[7:11],'linear','extrap');
ir_red=[power_red(1:6) in_val];
ir_blue=power_blue;
figure;scatter(1:11,power_blue,'Marker','+')
figure;scatter(1:11,ir_blue,'Marker','+')
hold on;scatter(1:11,ir_red,'Marker','+')
box off 


    
end