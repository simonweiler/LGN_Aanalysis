%Analysis of minis/failures and ramps for dLGN experiments (created by SW181211)

%%TR2019: this is just a test edit for pull request

%CHECK FLAG SECTION (line 14-21) AS WELL AS DIRECTORY SECTION (line 35-39)
%NESTED FUNCTIONS:  
%parseExperimentsXls_dLGN
%rampanalysis
%minianalyis,
%dLGN_plot_analysis

%%INITIATION

clear all;%delete all current variables in workspace
close all;%close all open windows/figures 

%%%%%%IMPORTANT FLAGS, PLEASE CHANGE HERE%%%%%%%%%%%%%%%%
analyze_mini=0;%flag if either mini only and/or ramp should be analyzed (1 or 0)
analyze_ramp=1;
fanalysis=0;
factor=4;%std threshold factor 
display=0;%flag to display plot (1 or 0)
ramp_rtrace=1;%save raw ephystraces or not (1 or 0)
savefile=1;%save file at the end or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analyze_mini==1 || analyze_ramp==1;
disp('dLGN Analysis Mini and Ramp');

sent = 'Which user data will be analyzed? type in 0 for SW or 1 for MF\n';%text appears in command window 
user = input(sent);%waiting for input which is either 0 or 1

experimentator = 'SW';%default SW data 
if user==1
experimentator = 'MF';%MF data 
end


dLGN_ephys={};%empty structure for saving variables
%%%%%%DIRECTORIES%%%%%%%
rdata_dir         = 'I:\Simon Weiler\EXPLORER ONE\dLGN_rawDATA';%data directory of raw data;change accordingly
adata_dir         = 'I:\Simon Weiler\EXPLORER ONE\dLGN_ephys_Analysis\';%data directory of extracted date;change accordingly 
ExpXls            = 'R:\Share\Simon\LGN_2019_SW_MF_JB_TR\dLGN_ephys_analysis_excel spread sheet\Experiments_dLGN_SW.xlsx';%directory where excel batch file is located;change accordingly 
%%%%%%%%%%%%%%%%%%%%%%%%

%% parse Experiments XLS database
batchopt          = parseExperimentsXls_dLGN(ExpXls,user);%calls the nested function parseExperimentsXls_dLGN and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed 
%% 


adder=1;%counting variable 
for i=1:nummice%for loop over experiments across days
  datapath=fullfile(rdata_dir, batchopt.mouse{i}, '\');%directory and name of experiments (from excel sheet)
  cd(char(datapath));%go to directory
 
  for k=1:length(batchopt.exp_ids{i})%loop in bigger loop for each cell per experimental day
      if batchopt.exp_ids{i}(k)<10%for cells with id less then XX0010, e.g., XX0001-XX0009
      n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));
      else
      n_str = sprintf( '%03d', batchopt.exp_ids{i}(k));%for cells with id mor then XX0010, e.g., XX0011-XX0099
      end
      fold_name=[experimentator n_str];%complete cell folder name such as SW0001 or MF0001
      exp_folder=fullfile(datapath,fold_name);%complete folder and directory
      list=dir([char(exp_folder) '\*.xsg']);%xsg files per cell 
      len=length(list);%number of xsg files per cell
      for j=1:len   
      load([char(exp_folder) '/' list(j).name],'-mat');%load each xsg file 
      iterations(:,j)=header.loopGui.loopGui.iterations;%find out whether mini or ramp recording 
      end
      
ramp=find(iterations==11);%ramp recordings 
failure1=find(iterations==50);%mini recordings 
failure2=find(iterations==100);%mini recordings 
disp(['CURRRENT EXPERIMENT is ', char(batchopt.mouse{i}), fold_name]);
if user==0
disp([num2str(length(ramp)/11),' ramp recordings']);
else
disp([num2str(length(ramp)),' ramp recordings']);
end
disp([num2str(length(failure1)),' failure recordings with 50 reps']);
disp([num2str(length(failure2)),' failure recordings with 100 reps']);
iterations=[];
%% RAMP ANALYSIS 
if analyze_ramp==1
  [blue_ramp, red_ramp]=rampanalysis(list, ramp, exp_folder, factor,display,ramp_rtrace,user);%use nested function rampanalysis 
  ramp=[];%clear variables for next iteration
  %list=[];%clear variables for next iteration
end

%% MINI ANALYSIS  
if analyze_mini==1
if length(failure1)>=1 & length(failure2)>=1%
[neg_failure, pos_failure PD1 PD2 IR1_r IR1_b IR2_b]=minianalysis(list, failure1, exp_folder, factor,display,user);%call minianalysis
[neg_failure, pos_failure PD1 PD2 IR1_r IR1_b IR2_b]=minianalysis(list, failure2, exp_folder, factor,display,user);%call minianalysis
failure1=[];
failure2=[];
%list=[];
elseif length(failure1)>=1 & length(failure2)==0%
 [neg_failure, pos_failure PD1 PD2 IR1_r IR1_b IR2_b]=minianalysis(list, failure1, exp_folder, factor,display,user);
failure1=[];
%list=[];
elseif length(failure2)>=1 & length(failure1)==0% 
[neg_failure, pos_failure PD1 PD2 IR1_r IR1_b IR2_b]=minianalysis(list, failure2, exp_folder, factor,display,user);
failure2=[];
else
disp('No failure recording');    
end
end

%prepare structure for all cells

if batchopt.exp_ids2{i}(k)==1
   slice_nr = 1 ;
else
   slice_nr = 2 ;
end
if analyze_ramp==1 & analyze_mini==1;  
 dLGN_ephys.data{1,1}='Animal ID';
 dLGN_ephys.data{1,2}='Experimental ID';
 dLGN_ephys.data{1,3}='Slice';
 dLGN_ephys.data{1,4}='Peak blue';
 dLGN_ephys.data{1,5}='Peak red';
 dLGN_ephys.data{1,6}='Category';
 dLGN_ephys.data{1,7}='Mini AMPA';
 dLGN_ephys.data{1,8}='Mini NMDA';
 dLGN_ephys.data{1,9}='Mini PD1';
 dLGN_ephys.data{1,10}='Mini PD2';
dLGN_ephys.data{1,11}='Mini Irradiance1';
 dLGN_ephys.data{1,12}='Mini Irradiance2';
 dLGN_ephys.data{1,13}='Mini blue const. Irradiance';
 dLGN_ephys.data{adder+1,1}=[char(batchopt.mouseID{i})];
 dLGN_ephys.data{adder+1,2}=[char(batchopt.mouse{i}), fold_name];
  dLGN_ephys.data{adder+1,3}=slice_nr;
 dLGN_ephys.data{adder+1,4}=blue_ramp;
 dLGN_ephys.data{adder+1,5}=red_ramp;
 dLGN_ephys.data{adder+1,6}=batchopt.exp_ids3{i}(k);
 
 dLGN_ephys.data{adder+1,7}=neg_failure;
 dLGN_ephys.data{adder+1,8}=pos_failure;
 dLGN_ephys.data{adder+1,9}=PD1;
 dLGN_ephys.data{adder+1,10}=PD2;
 dLGN_ephys.data{adder+1,11}=IR1_r;
 dLGN_ephys.data{adder+1,12}=IR1_b;
 dLGN_ephys.data{adder+1,13}=IR2_b;
 adder=adder+1;
end
if analyze_ramp==1 & analyze_mini==0;
 dLGN_ephys.data{1,1}='Animal ID';
 dLGN_ephys.data{1,2}='Experimental ID';
 dLGN_ephys.data{1,3}='Slice';
 dLGN_ephys.data{1,4}='Peak blue';
 dLGN_ephys.data{1,5}='Peak red';
 dLGN_ephys.data{1,6}='Category';
 dLGN_ephys.data{adder+1,1}=[char(batchopt.mouseID{i})];
 dLGN_ephys.data{adder+1,2}=[char(batchopt.mouse{i}), fold_name];
 dLGN_ephys.data{adder+1,3}=slice_nr;
 dLGN_ephys.data{adder+1,4}=blue_ramp;
 dLGN_ephys.data{adder+1,5}=red_ramp;
 dLGN_ephys.data{adder+1,6}=batchopt.exp_ids3{i}(k);
 adder=adder+1;
end
if analyze_ramp==0 & analyze_mini==1;
  dLGN_ephys.data{adder,1}=[char(batchopt.mouse{i}), fold_name];
  dLGN_ephys.data{adder,2}=neg_failure;
  dLGN_ephys.data{adder,3}=pos_failure;
  dLGN_ephys.data{adder,4}=PD1;
  dLGN_ephys.data{adder,5}=PD2;
  dLGN_ephys.data{adder,6}=IR1_r;
  dLGN_ephys.data{adder,7}=IR1_b;
  dLGN_ephys.data{adder,8}=IR2_b;
  dLGN_ephys.data{adder,9}=slice_nr;
  adder=adder+1;
end
end 
list=[];
end

% SAVE in analyzed directory   
if savefile==1
cd(adata_dir);
FileName=['Data_',experimentator,'_',datestr(now, 'hh-dd-mmm-yyyy')];
save(FileName,'-struct','dLGN_ephys');
disp('FILE SAVED');
else
disp('FILE NOT SAVED');
end
end
%% FURTHER ANALYSIS SUCH AS ODI or AMPA/NMDA RATIO

if fanalysis==1   
   disp('dLGN data plotting of extracted parameters and calculation of ODI and AMPA/NMDA RATIOS');
   adata_dir         = 'I:\Simon Weiler\EXPLORER ONE\dLGN_ephys_Analysis\';%data directory of saved data
   
   [ramps_peak ODI data]=dLGN_plot_analysis(adata_dir,0) 
   
end
  
    




    
   
    
%-------------------------------------------    

    




