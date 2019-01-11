function dLGN_plot_analysis(directory) 
%SW181229
%Function for claculating and plotting extracted dLGN ephys data. This
%function uses the extracted ephys data gnerated with script Analysis_mini_ramp as an input

%Inputs:
%directory: directory of extracted mat file gnerated with script Analysis_mini_ramp

%NESTED FUNCTION: 
%barwitherr
%uipickfiles

filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file 
%% RAMP ANALYSIS
%%LOAD PEAK AMPLITUDES FOR AMPA and NMDA 
for i=1:size(data,1)
%AMPA
blue_ramp_70(:,i)=data{i,2}.neg_peak2(1,:);
blue_laser_70(:,i)=data{i,2}.laser_amp(1,:);
red_ramp_70(:,i)=data{i,3}.neg_peak1(2,:);
red_laser_70(:,i)=data{i,3}.laser_amp(2,:); 
blue_constant_70(:,i)=data{i,2}.neg_peak2(2,:);
blue_claser_70(:,i)=data{i,2}.laser_amp(2,:);
%NMDA
blue_ramp_40(:,i)=data{i,2}.pos_peak2(3,:);
blue_laser_40(:,i)=data{i,2}.laser_amp(3,:);
red_ramp_40(:,i)=data{i,3}.pos_peak1(4,:);
red_laser_40(:,i)=data{i,3}.laser_amp(4,:);
blue_constant_40(:,i)=data{i,2}.pos_peak2(4,:);
blue_claser_40(:,i)=data{i,2}.laser_amp(4,:);
end
i=[];
%% 
%SET VALUES BELOW STD THRESHOLD CRITERION TO ZERO
for i=1:length(data)
%%%%%%%BLUE ONLY AMPA
idx=find(data{i,2}.neg_fail2(1,:)>0)';
d=data{i,2}.neg_peak2(1,:);
%%%%%%%
ze=zeros(11,1);
if length(idx)>=5;
ze(idx)=d(idx);
b_r_70(:,i)=ze ;
else
b_r_70(:,i)=ze;
end
%%%%%%%
idx=[];
d=[];
ze=[];
%%%%%%%RED AMPA
idx=find(data{i,3}.neg_fail1(2,:)>0)';
d=data{i,3}.neg_peak1(2,:);
%%%%%%%
ze=zeros(11,1);
if length(idx)>=5;
ze(idx)=d(idx);
r_r_70(:,i)=ze;
else
r_r_70(:,i)=ze;
end
%%%%%%%
idx=[];
d=[];
ze=[];
%%%%%%%BLUE CONSTANT AMPA
idx=find(data{i,2}.neg_fail2(2,:)>0)';
d=data{i,2}.neg_peak2(2,:);
%%%%%%%
ze=zeros(11,1);
if length(idx)>=5;
ze(idx)=d(idx);
b_c_70(:,i)=ze;
else
b_c_70(:,i)=ze;
end
%%%%%%%
idx=[];
d=[];
ze=[];
%%%%%%%BLUE ONLY NMDA
idx=find(data{i,2}.pos_fail2(3,:)>0)';
d=data{i,2}.pos_peak2(3,:);
%%%%%%%
ze=zeros(11,1);
if length(idx)>=5;
ze(idx)=d(idx);
b_r_40(:,i)=ze ;
else
b_r_40(:,i)=ze;
end
%%%%%%%
idx=[];
d=[];
ze=[];
%%%%%%%RED NMDA
idx=find(data{i,3}.pos_fail1(4,:)>0)';
d=data{i,3}.pos_peak1(4,:);
%%%%%%%
ze=zeros(11,1);
if length(idx)>=5;
ze(idx)=d(idx);
r_r_40(:,i)=ze;
else
r_r_40(:,i)=ze;
end
%%%%%%%
idx=[];
d=[];
ze=[];
%%%%%%%BLUE CONSTANT NMDA
idx=find(data{i,2}.pos_fail2(4,:)>0)';
d=data{i,2}.pos_peak2(4,:);
%%%%%%%
ze=zeros(11,1);
if length(idx)>=5;
ze(idx)=d(idx);
b_c_40(:,i)=ze;
else
b_c_40(:,i)=ze;
end
%%%%%%%
idx=[];
d=[];
ze=[];
end
%% 
%ODI and AMPA/NMDA Ratio
 for i=1:length(data)
 %AMPA
 ODI_raw_AMPA(:,i)=(abs(red_ramp_70(:,i))-abs(blue_constant_70(:,i)))./(abs(red_ramp_70(:,i))+abs(blue_constant_70(:,i)));
 ODI_AMPA(:,i)=(abs(r_r_70(:,i))-abs(b_c_70(:,i)))./(abs(r_r_70(:,i))+abs(b_c_70(:,i)));
 %NMDA
 ODI_raw_NMDA(:,i)=(abs(red_ramp_40(:,i))-abs(blue_constant_40(:,i)))./(abs(red_ramp_40(:,i))+abs(blue_constant_40(:,i)));
 ODI_NMDA(:,i)=(abs(r_r_40(:,i))-abs(b_c_40(:,i)))./(abs(r_r_40(:,i))+abs(b_c_40(:,i)));
 %AMPA/NMDA Ratio
 Ratio_an_r(:,i)=abs(r_r_70(:,i))./abs(r_r_40(:,i));
 Ratio_an_b(:,i)=abs(b_r_70(:,i))./abs(b_r_40(:,i));
 Ratio_an_bc(:,i)=abs(b_c_70(:,i))./abs(b_c_40(:,i));
 end
 %ODI average; used the last 6 entries from 11 ramps-> NEEDS TO BE DISCUSSED
 ODI_A=nanmean(ODI_AMPA(5:end,:));
 ODI_A_sem=nanstd(ODI_AMPA(5:end,:))/sqrt(length(data));
 ODI_N=nanmean(ODI_NMDA(5:end,:));
 ODI_N_sem=nanstd(ODI_NMDA(5:end,:))/sqrt(length(data));
 %RATIO AMPA/NMDA average 
 R_r=nanmean(Ratio_an_r(5:end,:));
 R_r_sem=nanstd(Ratio_an_r(5:end,:))/sqrt(length(data));
 R_b=nanmean(Ratio_an_b(5:end,:));
 R_b_sem=nanstd(Ratio_an_b(5:end,:))/sqrt(length(data));
 R_bc=nanmean(Ratio_an_bc(5:end,:));
 R_bc_sem=nanstd(Ratio_an_bc(5:end,:))/sqrt(length(data));
 %Create combined vector for plotting in bar graph later
 com_ODI=[ODI_A' ODI_N'];
 com_ODI_sem=[ODI_A_sem' ODI_N_sem'];
 com_R=[R_b' R_r' R_bc'];
 com_R_sem=[R_b_sem' R_r_sem' R_bc_sem'];
%% 
i=[];
%PLOT GRAPHS
cvec=unique(hsv(length(data)*30),'rows');
Legend=cell(length(data),1)%  two positions 
f1=figure('Name','AMPA');
set(gcf, 'Position', [200, 0, 600, 1500]);
f2=figure('Name','NMDA');
set(gcf, 'Position', [600, 0, 600, 1500]);

for i=1:length(data)
figure(f1);
subplot(3,1,1);
h{i}=plot(blue_laser_70(:,i),b_r_70(:,i),'--s','LineWidth',2,'MarkerSize',5,'Color',cvec(i*30,:));
axis square;
title('Blue laser only');
xlabel('Laser amplitude (a.u.)');
ylabel('Peak synaptic current (pA)');
Legend{i}=data{i,1};  
hold on;
end
figure(f1);
legend(Legend,'Location','EastOutside');
legend boxoff;
hold on;
i=[];
%%%second repeat
for i=1:length(data)
subplot(3,1,2);
m{i}=plot(blue_laser_70(:,i),r_r_70(:,i),'--s','LineWidth',2,'MarkerSize',5,'Color',cvec(i*30,:));
title('Red laser');
xlabel('Laser amplitude (a.u.)');
ylabel('Peak synaptic current (pA)');
axis square;
hold on;
end
figure(f1);
legend(Legend,'Location','EastOutside');
legend boxoff;
hold on;
i=[];
%%%third repeat
for i=1:length(data)
subplot(3,1,3);
m{i}=plot(blue_claser_70(:,i),b_c_70(:,i),'--s','LineWidth',2,'MarkerSize',5,'Color',cvec(i*30,:));
title('Blue constant');
xlabel('Laser amplitude (a.u.)');
ylabel('Peak synaptic current (pA)');
axis square;
hold on;
end
figure(f1);
legend(Legend,'Location','EastOutside');
legend boxoff;
hold on;
i=[];

%40mV NMDA
%1st repeat
figure(f2);
for i=1:length(data)
subplot(3,1,1);
m{i}=plot(blue_laser_40(:,i),b_r_40(:,i),'--s','LineWidth',2,'MarkerSize',5,'Color',cvec(i*30,:));
title('Blue laser only');
xlabel('Laser amplitude (a.u.)');
ylabel('Peak synaptic current (pA)');
axis square;
Legend{i}=data{i,1};
hold on;
end
figure(f2);
legend(Legend,'Location','EastOutside');
legend boxoff;
hold on;
i=[];
%2nd repeat
for i=1:length(data)
subplot(3,1,2);
m{i}=plot(blue_laser_40(:,i),r_r_40(:,i),'--s','LineWidth',2,'MarkerSize',5,'Color',cvec(i*30,:));
title('Red laser');
xlabel('Laser amplitude (a.u.)');
ylabel('Peak synaptic current (pA)');
axis square;
Legend{i}=data{i,1};
hold on;
end
figure(f2);
legend(Legend,'Location','EastOutside');
legend boxoff;
hold on;
i=[];
%3rd repeat
for i=1:length(data)
subplot(3,1,3);
m{i}=plot(blue_claser_40(:,i),b_c_40(:,i),'--s','LineWidth',2,'MarkerSize',5,'Color',cvec(i*30,:));
title('Blue constant');
xlabel('Laser amplitude (a.u.)');
ylabel('Peak synaptic current (pA)');
axis square;
Legend{i}=data{i,1};
hold on;
end
figure(f2);
legend(Legend,'Location','EastOutside');
legend boxoff;
hold on;
i=[];

%ODI AND RATIO AMPA/NMDA
%ODI
try
figure;barwitherr(com_ODI_sem,com_ODI);
catch
figure;bar(com_ODI);hold on;errorbar(com_ODI,com_ODI_sem,'.');
end
ylim([-1.3 1.3])
axis square;
xlabel('Cell');
ylabel('ODI');
title('Ocular dominance');
legend('AMPA','NMDA');
%RATIO AMPA/NMDA
try
figure;barwitherr(com_R_sem,com_R);
catch
figure;bar(com_R);hold on;errorbar(com_R,com_R_sem,'.');
end
axis square;
xlabel('Cell');
ylabel('Peak AMPA / Peak NMDA');
title('AMPA / NMDA Ratio');
legend('Blue only','Red','Blue constant');
%% MINI ANALYSIS
%AMPA
i_o_suc=zeros(1,size(data,1));
i_o_suc_sem=zeros(1,size(data,1));
for i=1:size(data,1)
mblue_70(:,i)=data{i,4}(1:50,1);
mred_70(:,i)=data{i,4}(1:50,2);
mblue_const_70(:,i)=data{i,5}(1:50,2);
PD_blue_70(:,i)=data{i,6}(1:50,1);
PD_red_70(:,i)=data{i,6}(1:50,2);
PD_blue_const_70(:,i)=data{i,7}(1:50,2);
PD_br_ratio(:,i)=mean(PD_blue_70(:,i))/mean(PD_red_70(:,i));
if length(find(mblue_70(:,i)<0))>5 & length(find(mred_70(:,i)<0))<7
    disp('IPSI ONLY');
    idx_success=find(mblue_70(1:50,i)<0);
    i_o_suc(:,i)=mean(mblue_70(idx_success));
    i_o_suc_sem(:,i)=std(mblue_70(idx_success,i))/length(sqrt(idx_success));    
else
    disp ('EITHER CONTRA ONLY OR BINOCULAR');
end

end



    
    
    

end
