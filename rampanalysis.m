function  [blue_ramp, red_ramp]=rampanalysis(list, idx, pathName, fc, show, ramp_rtrace, user);
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

%% TR2019: filtering
filterephys = 1;        % filtering yes/no?
cutoff      = 1000      % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4         % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Bessel'; % filter type ('Bessel' or 'Butter' (for Butterworth -> ). Default: Bessel. Use Bessel at > 4 order to prevent ripples)

if filterephys;
    disp('- - - - - - - -')
    disp(['Filtering: ' num2str(order) ' pole ' type '-Filter w/ ' num2str(cutoff) ' Hz cutoff']);
    disp('- - - - - - - -')
end

%create vector with start and end point for each ramp within the cell recording
%plot if wanted
if show==1
    figure;
    set(gcf, 'Position', [200, 0, 1500, 1000]);
end

if user==0%SW
    runramp=1:11:length(idx);
    runramp=[runramp runramp(end)+11];
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
            
            if filterephys % TR2019: filtering
                traces = lowpassfilt(traces, order, cutoff, sr, type);
            end
            
            photodiode=data.acquirer.trace_1;%photodiode (PD) signal
            try
                blue_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{4,counter+1}.amplitude;%blue laser amplitude set in ephus
            catch
                blue_amp(j,counter)=0;
            end
            try
                red_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{2,counter+1}.amplitude;%red laser amplitude set in ephus
            catch
                red_amp(j,counter)=0;
            end
            bs=traces(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
            bs_std=std(bs);%std of baseline trace
            bs_traces=traces-mean(traces(base_start*srF:base_end*srF,:));%subtract baseline
            bs_photodiode=photodiode-mean(photodiode(base_start*srF:base_end*srF,:));
            %for first window
            neg_peak1(j,counter)=min(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%negative peak within the red stimulation window
            pos_peak1(j,counter)=max(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%positive peak within the red stimulation window
            integ1(j,counter)=trapz(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%Integral within the red stimulation window
            neg_fail1(j,counter)=neg_peak1(j,counter)<fc*bs_std*(-1);%vector with binary values when neg peaks crossed definded std threshold
            pos_fail1(j,counter)=pos_peak1(j,counter)>fc*bs_std;%vector with binary values when pos peaks crossed definded std threshold
            %photodiode
            PD1(j,counter)=mean(bs_photodiode(redpeak_start*srF:redpeak_end*srF,:));%max values of PD signal within the red stimulation window
            %%%%extract irradiance for red%%%%
            yirr_red(j,counter)=(12.19*PD1(j,counter)-0.4319)/100;
            %for second window (same extraction as above for blue laser window
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%IMPLEMENTED AFTER MEETING FROM 190109%% neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
            %neg peak2 is calculated using the current difference  between the last 10ms of
            %the first time window and the peak in the subsequent 2nd window to
            %correct for decay issues from the first pulse
            neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:))-mean(bs_traces((redpeak_end-10)*srF:redpeak_end*srF,:));
            
            %%%IMPLEMENTED AFTER MEETING FROM 190109%For NMDA: approach is to fit an expontial and then subtract this from
            %the actual curve to detect a second peak
            if j<=2
                pos_peak2(j,counter)=0;
                pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                yf=bs_traces;
                diff_bs_traces=bs_traces;
            elseif j==3
                pos_peak2(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                yf=bs_traces;
                diff_bs_traces=bs_traces;
            else j==4;
                currmaxpos(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                if currmaxpos(j,counter)>pos_peak1(j,counter)
                    pos_peak2(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                    yf=bs_traces;
                    diff_bs_traces=bs_traces;
                    pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                else
                    xt=1:50000;
                    A=pos_peak1(j,counter);
                    t1=find(bs_traces==A);
                    t1=t1(1);
                    t=t1:redpeak_end*srF;
                    t=t';
                    curr_t=bs_traces(t);
                    try
                        [f gof]=fit(t,curr_t,'exp1');
                        yf=f.a*exp(f.b*xt);
                        for m=1:10000;
                            diff_bs_traces(m,:)=bs_traces(m)-yf(m);
                        end
                        bs_diff_std=std(diff_bs_traces((redpeak_end-100)*srF:redpeak_end*srF,:));
                        if show==1
                            %figure;plot(bs_traces);hold on;plot(yf);plot(diff_bs_traces);
                        end
                        pos_peak2(j,counter)=max(diff_bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                        pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_diff_std;
                        if gof.adjrsquare<0.9
                            pos_peak2(j,counter)=0;
                            pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_diff_std;
                        end
                    catch
                        pos_peak2(j,counter)=0;
                        pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                        yf=bs_traces;
                        diff_bs_traces=bs_traces;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            integ2(j,counter)=trapz(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
            neg_fail2(j,counter)=neg_peak2(j,counter)<fc*bs_std*(-1);
            %photodiode
            PD2(j,counter)=mean(bs_photodiode(bluepeak_start*srF:bluepeak_end*srF,:));
            %%%%extract irradiance for blue%%%%
            yirr_blue(j,counter)=(7.232*PD2(j,counter)-0.9951)/100;%given in mW/mm2 compare to Klapoetke 2014
            
            %ephys_traces
            ephys_traces(:,counter,j)=bs_traces;
            fit_traces(:,counter,j)=yf;
            diff_traces(:,counter,j)=diff_bs_traces;
            
            counter=counter+1;
            traces=[];
            %%%%%%%%%%%%%%plot
            if show==1
                plot(bs_traces(1:10000,:),'linewidth',1,'Color',[0 0 0]+0.05*counter);
                hold on;
                ylabel('Synaptic input (pA)');
                xlabel('Time (ms)');
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
    end
    %MF user==0
else
    for j=1:length(idx)% how many ramps in total; loop across ramps per cell
        counter=1;
        if show==1
            subplot(2,(length(idx))-2,j);
        end
        load([char(pathName) '/' list(idx(j)).name],'-mat');
        sr = header.ephys.ephys.sampleRate;%check sample rate
        srF = 1/(1000/sr);
        traces=data.ephys.trace_1;%raw ephys trace
        photodiode=data.acquirer.trace_1;%photodiode (PD) signal
        ind_traces=reshape(traces,[length(traces)/11 11]);
        photodiode=reshape(photodiode,[length(traces)/11 11]);
        for i=1:size(ind_traces,2);%within each ramp load xsg files (11 in total per ramp)
            traces_clip=ind_traces(:,i);
            photodiode_clip=photodiode(:,i);
            blue_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{2,counter+1}.amplitude;%blue laser amplitude set in ephus
            try
                red_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{3,counter+1}.amplitude;%red laser amplitude set in ephus
            catch
                red_amp(j,counter)=0;
            end
            bs=traces_clip(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
            bs_std=std(bs);%std of baseline trace
            bs_traces=traces_clip-mean(traces_clip(base_start*srF:base_end*srF,:));%subtract baseline
            bs_photodiode=photodiode_clip-mean(photodiode_clip(base_start*srF:base_end*srF,:));
            %for first window
            neg_peak1(j,counter)=min(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%negative peak within the red stimulation window
            pos_peak1(j,counter)=max(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%positive peak within the red stimulation window
            integ1(j,counter)=trapz(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%Integral within the red stimulation window
            neg_fail1(j,counter)=neg_peak1(j,counter)<fc*bs_std*(-1);%vector with binary values when neg peaks crossed definded std threshold
            pos_fail1(j,counter)=pos_peak1(j,counter)>fc*bs_std;%vector with binary values when pos peaks crossed definded std threshold
            %photodiode
            PD1(j,counter)=mean(bs_photodiode(redpeak_start*srF:redpeak_end*srF,:));%max values of PD signal within the red stimulation window
            %%%%extract irradiance for red%%%%
            yirr_red(j,counter)=(104.1 *PD1(j,counter)-3.467)/100;
            %for second window (same extraction as above for blue laser window
            %     neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
            %     pos_peak2(j,counter)=max(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%IMPLEMENTED AFTER MEETING FROM 190109%% neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
            %neg peak2 is calculated using the current difference  between the last 10ms of
            %the first time window and the peak in the subsequent 2nd window to
            %correct for decay issues from the first pulse
            neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:))-mean(bs_traces((redpeak_end-10)*srF:redpeak_end*srF,:));
            
            %%%IMPLEMENTED AFTER MEETING FROM 190109%For NMDA: approach is to fit an expontial and then subtract this from
            %the actual curve to detect a second peak
            if j<=2
                pos_peak2(j,counter)=0;
                pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                yf=bs_traces;
                diff_bs_traces=bs_traces;
            elseif j==3
                pos_peak2(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                yf=bs_traces;
                diff_bs_traces=bs_traces;
            else j==4;
                currmaxpos(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                if currmaxpos(j,counter)>pos_peak1(j,counter)
                    pos_peak2(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                    yf=bs_traces;
                    diff_bs_traces=bs_traces;
                    pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                else
                    xt=1:200000;
                    A=pos_peak1(j,counter);
                    t1=find(bs_traces==A);
                    t1=t1(1);
                    t=t1:redpeak_end*srF;
                    t=t';
                    curr_t=bs_traces(t);
                    try
                        [f gof]=fit(t,curr_t,'exp1');
                        yf=f.a*exp(f.b*xt);
                        for m=1:40000;
                            diff_bs_traces(m,:)=bs_traces(m)-yf(m);
                        end
                        bs_diff_std=std(diff_bs_traces((redpeak_end-100)*srF:redpeak_end*srF,:));
                        if show==1
                            %figure;plot(bs_traces);hold on;plot(yf);plot(diff_bs_traces);
                        end
                        pos_peak2(j,counter)=max(diff_bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
                        pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_diff_std;
                        if gof.adjrsquare<0.9
                            pos_peak2(j,counter)=0;
                            pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_diff_std;
                        end
                    catch
                        pos_peak2(j,counter)=0;
                        pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
                        yf=bs_traces;
                        diff_bs_traces=bs_traces;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            integ2(j,counter)=trapz(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
            neg_fail2(j,counter)=neg_peak2(j,counter)<fc*bs_std*(-1);
            % pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std;
            %photodiode
            PD2(j,counter)=mean(bs_photodiode(bluepeak_start*srF:bluepeak_end*srF,:));
            %%%%extract irradiance for blue%%%%
            yirr_blue(j,counter)=(679.2*PD2(j,counter)-26.82)/100;
            %ephys_traces
            ephys_traces(:,counter,j)=bs_traces;
            fit_traces(:,counter,j)=yf;
            diff_traces(:,counter,j)=diff_bs_traces;
            
            counter=counter+1;
            traces=[];
            
            %%%%%%%%%%%%%%plot
            if show==1
                plot(bs_traces(1:20000,:),'linewidth',1,'Color',[0 0 0]+0.05*counter);
                hold on;
                ylabel('Synaptic input (pA)');
                xlabel('Time (ms)');
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
    end
end
%%%%%%%%%%%%%%%%%%%%%  output %%%%%%%%%%%%%%%%%
red_ramp.neg_peak1=neg_peak1;
red_ramp.pos_peak1=pos_peak1;
red_ramp.integ1=integ1;
red_ramp.neg_fail1=neg_fail1;
red_ramp.pos_fail1=pos_fail1;
red_ramp.PD=PD1;
red_ramp.irr_red=yirr_red;
red_ramp.laser_amp=red_amp;
if ramp_rtrace==1;
    red_ramp.ephys_traces=ephys_traces;
    red_ramp.fit_traces=fit_traces;
    red_ramp.diff_traces=diff_traces;
    
end

%create structure with extracted parameters
blue_ramp.neg_peak2=neg_peak2;
blue_ramp.pos_peak2=pos_peak2;
blue_ramp.integ2=integ2;
blue_ramp.neg_fail2=neg_fail2;
blue_ramp.pos_fail2=pos_fail2;
blue_ramp.PD=PD2;
blue_ramp.irr_blue=yirr_blue;
blue_ramp.laser_amp=blue_amp;

end
