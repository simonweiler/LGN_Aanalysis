function r=morpho_LGN(morpho_folder,morpho_data,xf,yf,zf)

%% LOAD

len_br=length(morpho_data);
cd(morpho_folder);
 for m=1:len_br
dendr(m)=load_tree(morpho_data(m).name);
dendr(m).X=dendr(m).X/(xf(i));
dendr(m).Y=dendr(m).Y/(yf(i));
dendr(m).Z=dendr(m).Z/(zf(i));
X{m}=dendr(m).X;
Y{m}=dendr(m).Y;
Z{m}=dendr(m).Z;
% hold on;
% scatter3(dendr(m).X,dendr(m).Y,dendr(m).Z,'.');
% hold on;
names{m}=dendr(m).name; 
 end
 Index = find(contains(names,'Soma'));
 %% CONCATENATE 

structure_map = zeros(size(dendr));
%for all the cells
    for k=1:size(dendr,2)-1
        %check whether the field is empty
        if ~isempty(dendr(k).dA)
            %label the field in the map
            structure_map(k) = 1;
            %save the position
            nonempty_pos = k;
        end
    %label the last non-empty pos
    structure_map(nonempty_pos) = 2;
    end
    cat_fields = dendr(structure_map(1,:)==2)
    %get the number of fields
    num_fields = sum(structure_map(1,:)==2);
    %for all the involved fields
    for fields = 2:num_fields
        %if it's the first iteration
        if fields == 2
            %concatenate the first 2 trees
            result_tree = cat_tree(cat_fields(1),cat_fields(2));
        else
            %keep concatenating the fields
            result_tree = cat_tree(result_tree,cat_fields(fields));
        end
        %transpose the R vector in result_tree
        result_tree.R = result_tree.R';
    end
  
    %rename the name field
    result_tree.name = 'All';
    
    %get the soma structure and subtract it from the rest
    soma_struct = dendr(1,structure_map(1,:)==0); 
    mx=mean(soma_struct.X);
    my=mean(soma_struct.Y);
    mz=mean(soma_struct.Z);
    %create tree with soma subtracted
    result_tree.X= result_tree.X-mx;
    result_tree.Y= result_tree.Y-my;
    result_tree.Z= result_tree.Z-mz;
   somasub=soma_struct;
   somasub.X=somasub.X-mx;
   somasub.Y=somasub.Y-my;
   somasub.Z=somasub.Z-mz;
   
   %PLOT
    figure;
    scatter3(result_tree.X,result_tree.Y,result_tree.Z,'.');
    hold on;
	somasub.X,somasub.Y,somasub.Z,'.');
    
    traces=result_tree;
    %save the resulting tree in a cell
stats=stats_tree(traces,'-w -x');%TREES toolbox, download functions and put in Matlab path! check and use search find in the manual : https://www.treestoolbox.org/downloads/TREES_manual.pdf
traX=traces.X;%x coordinates
traY=traces.Y;%y coordinates
traZ=traces.Z;%z coordinates
euc=eucl_tree(traces); %euclidean distances of nodes to root [um]
len=len_tree(traces);% vector containing length values of tree segments [um]
path_len=Pvec_tree (traces, len);% path length from the root (µm)
tr=repair_tree(traces);%get rid of trifurcations 
    
%%%%%%%%%%%%%
%check ABI Supplementary Table 3 list : https://www.biorxiv.org/content/early/2018/07/17/368456

%Branching Pattern Feature;  
r.depth=abs(min(traces.Z)-max(traces.Z));%P2: single value; max Z subtracted from min Z
r.early_branch=min(length(stats.dstats.blen))/stats.gstats.max_plen;%P3: single value; not clear what it means, SW divided minimum branch length by max path length
r.height=abs(min(traces.Y)-max(traces.Y));%P4: single value; max Y subtracted from min Y
r.high_xyz=max([traces.X traces.Y traces.Z]);%P5: vector; not sure
r.low_xyz=min([traces.X traces.Y traces.Z]);%P6: vector; not sure 
r.max_branch_order=stats.gstats.maxbo;%P7: single value 
r.max_euclidean_distance=max(euc);%P8: single value 
r.max_path_distance=stats.gstats.max_plen;%P9: single value    
%P10, a little bit of precalculations
idx=find(typeN_tree(traces)==0 | typeN_tree(traces)==2);%indices for bifurcation (2) or terminal point (0)
summed_euc_bif=sum(abs(diff(euc(idx))));%sum of the euclidean difference between bifurcations and/or tips
summed_path_bif=sum(abs(diff(path_len(idx))));%sum of the path distance difference between bifurcations and/or tips
r.mean_contraction=summed_euc_bif/summed_path_bif;%P10; ratio 
r.num_bifurcations=length(find(B_tree(traces)));%P13 single value
r.num_branches=length(stats.dstats.blen);%P14 single value     
r.num_stems=len_br-1;
r.num_tips=sum(T_tree(traces));%P19 single value; number of terminal points 
r.total_length=stats.gstats.len;%P21 single value 
r.width=abs(min(traces.X)-max(traces.X));%P24: single value; max X subtracted from min X
r.bifurcation_angle_local=nanmean(abs(angleB_tree (tr)));%P25: single value
r.bifurcation_kurt_xyz=[kurtosis(traX(find(B_tree(traces)))) kurtosis(traY(find(B_tree(traces)))) kurtosis(traZ(find(B_tree(traces))))];%P27: vector with three entries for kurtosis of x,y and z
r.bifurcation_skew_xyz=[skewness(traX(find(B_tree(traces)))) skewness(traY(find(B_tree(traces)))) skewness(traZ(find(B_tree(traces))))];%P28: vector with three entries for skewness of x,y and z
r.bifurcation_std_xyz=[std(traX(find(B_tree(traces)))) std(traY(find(B_tree(traces)))) std(traZ(find(B_tree(traces))))];%P29: vector
%P30
com=[traX(find(B_tree(traces))) traY(find(B_tree(traces))) traZ(find(B_tree(traces)))];
r.first_bifurcation_moment_xyz=table2array(regionprops3(true(size(com)), com, 'WeightedCentroid'));%P30; not sure
r.second_bifurcation_moment_xyz=[var(traX(find(B_tree(traces)))) var(traY(find(B_tree(traces)))) var(traZ(find(B_tree(traces))))];%P31; not sure
    
%Combined Features 
r.density=r.total_length/(r.widthr.*height);%P40 single value ; not sure
r.bifurcation_centroid_over_distance_xyz=r.first_bifurcation_moment_xyz./[r.width r.height r.depth]%P42 vector 
r.bifurcation_centroid_over_stdev=r.first_bifurcation_moment_xyz./sqrt(r.second_bifurcation_moment_xyz);%P43 vector
r.bifurcation_stdev_over_centroid_xyz=r.bifurcation_std_xyz./r.first_bifurcation_moment_xyz;%P44 vector
r.bifurcation_stdev_over_distance_xyz=r.bifurcation_std_xyz./[r.width r.height r.depth];%P45 vector     
    
%Sholl Analysis 
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(traces, 20, '-s');

a1x=[0:150]*-1;
a1y=[0:150];

b1x=[0:150];
b1y=[0:150];

a2x=[0:150];
a2y=[0:150]*-1;

b2x=[0:150]*-1;
b2y=[0:150]*-1;
% 
% a1_idx=(XP >= a1x(end)) & (XP <= a1x(1)) & (YP >= a1y(1)) & (YP <= a1y(end))
% a2_idx=(XP >= a2x(1)) & (XP <= a2x(end)) & (YP >= a2y(end)) & (YP <= a2y(1))
% b1_idx=(XP >= b1x(1)) & (XP <= b1x(end)) & (YP >= b1y(1)) & (YP <= b1y(end))
% b2_idx=(XP >= b2x(end)) & (XP <= b2x(1)) & (YP >= b2y(end)) & (YP <= b2y(1))
a1_counts= sum((XP >= a1x(end)) & (XP <= a1x(1)) & (YP >= a1y(1)) & (YP <= a1y(end)));
a2_counts= sum((XP >= a2x(1)) & (XP <= a2x(end)) & (YP >= a2y(end)) & (YP <= a2y(1)));
b1_counts= sum((XP >= b1x(1)) & (XP <= b1x(end)) & (YP >= b1y(1)) & (YP <= b1y(end)));
b2_counts= sum((XP >= b2x(end)) & (XP <= b2x(1)) & (YP >= b2y(end)) & (YP <= b2y(1)));

com_a=a1_counts+ a2_counts;
com_b=b1_counts+ b2_counts;
com_planes=[com_a com_b];
r.DOi=min(com_planes)/max(com_planes);% from Krahe et al. 2011 JNeuroscience

%VISUALIZATION OF DENDRTITIC FIELD
[theta,rho] = cart2pol(XP,YP);
% for g=1:length(dd);
%   jj=find(iD==g) 
%   if jj>0
%  rho_outer=max(rho(jj));
%  rho_int(g)=max(rho(jj));
%   elsefff
%      rho_outer=0;
%      rho_int(g)=0;
%   end
%   mh=find(rho==rho_int(g));
%   if ~isempty(mh)
%   rho_idx=mh(end);
%   else 
%   rho_idx=1;
%   end
%   theta_int(g)=theta(rho_idx);
%   
% end
figure;h = polarhistogram(theta,20)
h.DisplayStyle = 'stairs';

  
r.com=[r.apical_mradial r.apical_len r.apical_plen r.apical_bpoints r.apical_maxbo r.apical_mblen r.apical_mplen r.apical_wh r.apical_xspan r.apical_yspan r.basal_mradial r.basal_len r.basal_plen r.basal_bpoints r.basal_maxbo r.basal_mblen r.basal_mplen r.basal_wh r.basal_xspan r.basal_yspan r.leng];
r.traces=[r.apical_trees;r.basal_trees;r.soma;r.cellnames]';
end


