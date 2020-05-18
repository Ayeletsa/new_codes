%% Population vectore analysis
clear

%% data:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_info=dir(dir_data);
analysis_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\PV_analysis\';
%% load parameters
population_vector_param_file_name='D:\Ayelet\2bat_proj\Analysis\new_code\params\population_vector_params.mat';
load(population_vector_param_file_name);
pop_vec_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\pop_vec\';

prm.fields.bin_size=allo_X_bin_size;
prm.fields.bin_limits=allo_X_bins_vector;
prm.fields.min_time_spent_per_meter=allo_time_spent_minimum_for_1D_bins;
prm.fields.ker_SD=ker_SD;
prm.fields.pos_fs=frames_per_second;
%% initialize:
cell_dir_count=0;

PV_co_dir_bin=[];
PV_solo_dir=[];
day=0;
day_i=0;
file_names={dir_info.name};
file_names=file_names(find([dir_info.isdir]==0));

%% 1. run over cells:
for cell_i=1:length(file_names)
%for cell_i=1:30
    % load data:
    load(fullfile(dir_data,file_names{cell_i}))
    % go over conditions to run the analysis:
    for dir_i=1:2
        a=sum(~isnan(cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(dir_i).spikes.ts_usec(:)))>=min_n_spike;
        b=cell_co_solo_initial_analysis.solo(dir_i).SI>SI_threshold;
        dir_cond(dir_i)=a & b;
    end
    
    if dir_cond(1) | dir_cond(2)
        %%  load data of both dirs to struct:
        %1. read cell data
        data=PV_get_cell_data(cell_co_solo_initial_analysis);
        %2. run on intersect data of y:
        if run_only_on_intersent_data==1
            data=PV_get_y_intersect_data(data,allo_X_bins_vector,allo_X_bins_vector_of_centers);
        end
        
        %% run over directions:
        for dir_i=1:2
            
            if dir_cond(dir_i)
                cell_dir_count=cell_dir_count+1;
                prev_day=day;
                day=cell_co_solo_initial_analysis.exp_data.day;
                
                %% compute allocentric tuning curve during solo in both directions
                solo_bsp=data(dir_i).solo.bsp.flight_x_pos(:);
                solo_spikes=data(dir_i).solo.spikes.flight_x_pos(:);
                solo_bsp=solo_bsp(~isnan(solo_bsp));
                solo_spikes=solo_spikes(~isnan(solo_spikes));
                [PSTH_solo(:,cell_dir_count),SI_solo(:,cell_dir_count)]=PV_arange_data_and_compute_shuffle_PSTH(solo_bsp,solo_spikes,prm);
                %other dir solo psth
                solo_bsp=data(3-dir_i).solo.bsp.flight_x_pos(:);
                solo_spikes=data(3-dir_i).solo.spikes.flight_x_pos(:);
                solo_bsp=solo_bsp(~isnan(solo_bsp));
                solo_spikes=solo_spikes(~isnan(solo_spikes));
                [PSTH_solo_other_dir(:,cell_dir_count),SI_solo_other_dir(:,cell_dir_count)]=PV_arange_data_and_compute_shuffle_PSTH(solo_bsp,solo_spikes,prm);
              
                %% comupte tuning curves for co data bin by bin of CO:
                n_co=size(data(dir_i).co.bsp.flight_ts,1);
                for bin_i=1:length(ego_bin_start)
                    bin_i
                    
                    %1. compute PSTH of CO data per egocentric bin:
                    %-----------------------------------------------
                    bsp_data=data(dir_i).co.bsp;
                    spike_data=data(dir_i).co.spikes;
                    
                    [FE,data_for_shuffle]=PV_find_relevant_data_per_bin(bsp_data,spike_data,ego_bin_start,ego_bin_end,bin_i,n_co);
                    FE_PSTH = FE_compute_PSTH(FE,prm);
                    PSTH_co_ego_bins(:,cell_dir_count,bin_i)=FE_PSTH.PSTH;
                   
                    SI_co_ego_bins(cell_dir_count,bin_i)=FE_PSTH.SI_bits_spike;
                    
                    %2. compute shuffle tuning curve
                    %-----------------------------------------------
                    %shuffle here will be solo tuning curve with data with
                    %the same range as in co:
                    
                    parfor shuffle_i=1:n_shuffles
                        shuffle_i
                        % a. find data for shuffle:
                        solo_data=data(dir_i).solo;
                        co_data=data(dir_i).co;
                        [shuffle_solo_bsp,shuffle_solo_spikes]=PV_find_solo_data_for_shuffle(solo_data,co_data,ego_bin_start,ego_bin_end,bin_i,dir_i) ;
                        shuffle_solo_bsp=shuffle_solo_bsp(~isnan(shuffle_solo_bsp));
                        shuffle_solo_spikes=shuffle_solo_spikes(~isnan(shuffle_solo_spikes));
                        % compute tuning curve:
                        [PSTH_solo_shuffle_ego_bins(:,cell_dir_count,shuffle_i,bin_i),SI_solo_shuffle_ego_bins(cell_dir_count,shuffle_i,bin_i)]=PV_arange_data_and_compute_shuffle_PSTH(shuffle_solo_bsp,shuffle_solo_spikes,prm);
                    end
                    
                end
                
                %% analyze behavior:
                if day~=prev_day
                    day_i=day_i+1;
                   [intersect_by_union(:,day_i),y_dev_per_pos(:,day_i),vel(:,day_i)]=PV_compute_behav(data,allo_X_bins_vector,ego_bin_start,ego_bin_end,dir_i);
                 
                end
                
            end
        end
        
    end
    
    
end
%% compute PV corr:
corrs=PV_compute_all_corrs(ego_bin_start,PSTH_co_ego_bins,PSTH_solo,PSTH_solo_shuffle_ego_bins,PSTH_solo_other_dir,n_shuffles);  

%% save data:
PV_data.corrs=corrs;
PV_data.PSTH_co_ego_bins=PSTH_co_ego_bins;
PV_data.solo_PSTH=PSTH_solo;
PV_data.PSTH_solo_shuffle_ego_bins=PSTH_solo_shuffle_ego_bins;
PV_data.solo_PSTH_other_dir=PSTH_solo_other_dir;

file_name=fullfile(analysis_struct_folder,'PV_mats_and_data.mat');
save(file_name, '-struct', 'PV_data')


%% plot:
%=============================================================================
figure('units','normalized','outerposition',[0 0 1 1])
%pos:
width=.1;
height=.1;
hor_dis=0.2;
ver_dis=0.18;

x_pos(1)=0.08;
x_pos(2)=x_pos(1)+hor_dis;
x_pos(3)=x_pos(2)+hor_dis;
x_pos(4)=x_pos(3)+hor_dis;
x_pos(5)=x_pos(4)+hor_dis;
x_pos(6)=x_pos(5)+hor_dis;
x_pos(7)=x_pos(6)+hor_dis;
x_pos(8)=x_pos(7)+hor_dis;

y_pos(1)=0.8;
y_pos(2)=y_pos(1)-ver_dis*0.8;
y_pos(3)=y_pos(2)-ver_dis;
y_pos(4)=y_pos(3)-ver_dis*0.8;
y_pos(5)=y_pos(4)-ver_dis;
y_pos(6)=y_pos(5)-ver_dis;

%1. plot cell cor:
%=======================================
pos_tuning=[x_pos(1),y_pos(1),width,height];
pos_mat=[x_pos(1),y_pos(2),width,height];
relevant_corr=corrs.cells;
y_label='cells';
title_str='corr of PSTH per cell averaged across cells';

plot_corr_mat_and_tuning(relevant_corr,pos_tuning,pos_mat,ego_bin_center,ego_bin_dis,y_label,title_str)

%2. plot PV corr:
%=======================================
pos_tuning=[x_pos(2),y_pos(1),width,height];
pos_mat=[x_pos(2),y_pos(2),width,height];
relevant_corr=corrs.PV;
y_label='x pos in tunnel';
title_str='Corr of allo population vector\nper position averaged across positions\n- no normalization';
plot_corr_mat_and_tuning(relevant_corr,pos_tuning,pos_mat,ego_bin_center,ego_bin_dis,y_label,title_str)

%3. plot PV corr norm to max:
%=======================================
pos_tuning=[x_pos(3),y_pos(1),width,height];
pos_mat=[x_pos(3),y_pos(2),width,height];
relevant_corr=corrs.norm_PV;
y_label='x pos in tunnel';
title_str='corr of population vector per position - norm to max';

plot_corr_mat_and_tuning(relevant_corr,pos_tuning,pos_mat,ego_bin_center,ego_bin_dis,y_label,title_str)

%4. plot SI
%=======================================
axes('position',[x_pos(4),y_pos(1),width,height]);hold on;
shuffle_SI_mean=squeeze(nanmean(SI_solo_shuffle_ego_bins,1));
co_SI_mean=nanmean(SI_co_ego_bins,1);
plot(ego_bin_center,shuffle_SI_mean,'color',[.5 .5 .5],'linewidth',2)
hold on;
plot(ego_bin_center,co_SI_mean,'linewidth',2,'color','m')
ylim([0 max([max(co_SI_mean)';max(shuffle_SI_mean(:))'])])
ylabel('SI')
xlabel('Inter-bat distance (m)')

%5. plot velocity
%=======================================
axes('position',[x_pos(1),y_pos(5),width,height])
plot(ego_bin_center,nanmean(vel,2),'linewidth',2,'color','m')
ylim([0 max(nanmean(vel))*1.05])
xlim([min(ego_bin_dis) max(ego_bin_dis)])
xlabel('Inter-bat distance (m)')
ylabel('Velocity (m/s)')

%5. plot deviation in y
%=======================================
axes('position',[x_pos(2),y_pos(5),width,height])
plot(ego_bin_center,nanmean(y_dev_per_pos,2),'linewidth',2,'color','m')
ylim([0 max(nanmean(y_dev_per_pos,2))*1.05])
xlim([min(ego_bin_dis) max(ego_bin_dis)])
xlabel('Inter-bat distance (m)')
ylabel('Deviation in y')

%6. Intersection/union
%=======================================
intersect_by_union(intersect_by_union==0)=nan;
axes('position',[x_pos(3),y_pos(5),width,height])
[h x]=hist(intersect_by_union)
bar(x,h/sum(h))
title(['full data intersect/union',num2str(nanmean(intersect_by_union(:)))])
xlim([0 1])

%save figure
fig_name=fullfile(pop_vec_folder,['pop_vector_fig_dir_',num2str(dir_i),'_binsize_',num2str(allo_X_bin_size),'min_spike',num2str(min_n_spike),'intersect_data_',num2str(run_only_on_intersent_data),'.png']);
saveas(gcf,fig_name)
%%
function plot_corr_mat_and_tuning(relevant_corr,pos_tuning,pos_mat,ego_bin_center,ego_bin_dis,y_label,title_str)
axes('position',pos_tuning);hold on;
plot(ego_bin_center,relevant_corr.co.shuffle.tuning,'color',[.5 .5 .5])
plot(ego_bin_center,relevant_corr.co.tuning,'linewidth',2,'color','m')
plot(ego_bin_center,relevant_corr.co.lower_bound.tuning,'color',[.5 .5 .5])
plot(ego_bin_center,relevant_corr.other_dir.tuning,'color',[0 0 0])

ylabel('r')
xlabel('Inter-bat distance (m)')
title(sprintf(title_str))

axes('position',pos_mat)
normalize_to_other_max=[];
cells_bin=1:size(relevant_corr.co.mat,1);
cells_bin_vector_of_centers=0.5:0.5:size(relevant_corr.co.mat,1)-0.5;
[xlimits, ylimits] = fn_plot_2D_field (relevant_corr.co.mat, ego_bin_dis, ego_bin_center,cells_bin, cells_bin_vector_of_centers,normalize_to_other_max);
xlabel('Inter-bat distance (m)')
ylabel(y_label)

end

