
%function to plot somt of the solo params
%1. check different kernel sizes:
k_sizes=[0.5,1,1.5];
%% figure prop:
figure('units','normalized','outerposition',[0 0 1 1])
x_pos=[0.05 0.55];
y_pos=[0.05 0.35 0.65];
raster_size=[0.2 0.15];
PSTH_size=[0.2 0.08];
fsize=10;
%%
figure_output_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\test_solo_params\';
data_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis_k_0.5\';
general_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis_k_';
dir_info=dir(data_folder_name);
behavior_struct_names = {dir_info.name};

SI_threshold=1;
%% params for valid cells during CO
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
param_file_name=fullfile(param_folder,'co_population_params.mat');
load(param_file_name)

for cell_i=3:length(behavior_struct_names)
    %load data of first kernel
    %% load data
    
    struct_name =behavior_struct_names{cell_i};
    file_name = fullfile(data_folder_name,struct_name);
    load(file_name);
    behavior_struct=cell_co_solo_initial_analysis;
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    cell_num=behavior_struct.exp_data.cell_num;
    
    % check if include cell:
    a=max([behavior_struct.solo.SI])>SI_threshold;
    b=(sum(~isnan(cell_co_solo_initial_analysis.solo(1).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(1).spikes.ts_usec(:))))>min_n_spike;
    c=(sum(~isnan(cell_co_solo_initial_analysis.solo(2).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(2).spikes.ts_usec(:))))>min_n_spike;
    d=b | c;
    e=cell_co_solo_initial_analysis.exp_data.mean_fr<max_for_pyramidal;
    cond_vec=[a,d,e];
    if sum(cond_vec)==length(cond_vec)
    
    %figure title
    % write cell's ID
    str = sprintf('Cell #%d, Date: %d, Bat %d',cell_num,day,bat) ;
    annotation('textbox',[.06 .88 .3 .1],'string',str,'EdgeColor','none','fontsize',15)
    for ii=1:length(k_sizes)
            for ii_dir=1:2

        % load data
        if ii~=1
            k_folder_name=[general_folder_name,num2str(k_sizes(ii)),'\'];
            file_name = [k_folder_name,'\bat',num2str(bat),'_day_',num2str(day),'_cell_', num2str(cell_num),'.mat'];
            load(file_name);
            behavior_struct=cell_co_solo_initial_analysis;
            
        end
        raster_pos=behavior_struct.solo(ii_dir).spikes.x_pos(:);
        raster_ts=behavior_struct.solo(ii_dir).spikes.ts_usec(:);
        time_limits = [min(behavior_struct.solo(ii_dir).bsp.ts_usec(:)) max(behavior_struct.solo(ii_dir).bsp.ts_usec(:))];
        bins=behavior_struct.solo(1).x_pos_firing_rate{1, 2};
        PSTH=behavior_struct.solo(ii_dir).PSTH_for_field_detection;
        fields_height=behavior_struct.solo(ii_dir).field_height;
        fields_center=behavior_struct.solo(ii_dir).field_center;
        plot_pos_1=[x_pos(ii_dir), y_pos(ii),raster_size];
        plot_pos_2=[x_pos(ii_dir), y_pos(ii)+raster_size(2),PSTH_size];
        title_str=sprintf('kernel=%.1f',k_sizes(ii));
        plot_solo_rater(raster_pos,raster_ts,time_limits,bins,PSTH,fields_height,fields_center,plot_pos_1,plot_pos_2,title_str)
    end
    end
    % save figure
    figure_name = [figure_output_folder,'\solo_params_bat',num2str(bat),'_day_',num2str(day),'_cell_', num2str(cell_num),'.jpg'];
    saveas(gcf,figure_name)
    clf
    end
end
%%
function plot_solo_rater(raster_pos,raster_ts,time_limits,bins,PSTH,fields_height,fields_center,plot_pos_1,plot_pos_2,title_str)
tunnel_limits=[0 135];
% raster of solo spikes
axes('position',plot_pos_1)
hold on
plot(raster_pos,raster_ts,'.','color','k','markersize',8)
time_labels=round(round(time_limits*1e-6*(1/60)*10-time_limits(1)*1e-6*(1/60)*10)/10);
set(gca,'ylim',time_limits,'ytick',time_limits,'yticklabel',time_labels,'xlim',tunnel_limits,'xtick',tunnel_limits)
ylabel ({'Time';'(min)'})
xlabel('pos (m)')

% solo firing rate - regular bins
axes('position',plot_pos_2)
hold on
plot(bins ,PSTH  ,'color','k','LineWidth',2);
plot(fields_center,fields_height,'*r')
max_x = 1.1*ceil(max(PSTH));
if max_x~=0 & ~isnan(max_x)
    set(gca,'xlim',tunnel_limits,'XTick',[],'ylim',[0 max_x],'ytick',max_x);
end
str = sprintf('Hz');
ylabel(str)
box off
title(title_str)
end

