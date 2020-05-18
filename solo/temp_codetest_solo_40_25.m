
%function to plot somt of the solo params
%1. check different kernel sizes:
k_sizes=[0.5,1,1.5];
th_vec=[0.5,1];
%% figure prop:
figure('units','normalized','outerposition',[0 0 1 1])
x_pos=[0.1 0.5];
y_pos=[0.35 0.65];
raster_size=[0.19 0.15];
PSTH_size=[0.19 0.08];
fsize=10;
dir_pos=[0, 0.5];
%%
figure_output_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\test_solo_params\';
data_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\co_40_solo_25\';
solo_40_foder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
%general_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis_k_';
dir_info=dir(data_folder_name);
behavior_struct_names = {dir_info.name};
behavior_struct_names=behavior_struct_names(find([dir_info.isdir]==0));
SI_threshold=1;
min_n_spike=100;
%% params for valid cells during CO
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
param_file_name=fullfile(param_folder,'co_population_params.mat');
load(param_file_name)

for cell_i=1:length(behavior_struct_names)
    %load data of first kernel
    %% load data
    
    struct_name =behavior_struct_names{cell_i};
    file_name = fullfile(data_folder_name,struct_name);
    load(file_name);
    solo_25=cell_co_solo_initial_analysis;
    file_name = fullfile(solo_40_foder_name,struct_name);
    load(file_name);
    solo_40=cell_co_solo_initial_analysis;
    
    bat=solo_25.exp_data.bat;
    day=solo_25.exp_data.day;
    cell_num=solo_25.exp_data.cell_num;
    
    % check if include cell:
    a=max([solo_25.solo.SI])>SI_threshold;
    b=(sum(~isnan(solo_25.solo(1).spikes.ts_usec(:)))+sum(~isnan(solo_25.co(1).spikes.ts_usec(:))))>min_n_spike;
    c=(sum(~isnan(solo_25.solo(2).spikes.ts_usec(:)))+sum(~isnan(solo_25.co(2).spikes.ts_usec(:))))>min_n_spike;
    d=b | c;
    e=solo_25.exp_data.mean_fr<max_for_pyramidal;
    cond_vec=[a,d,e];
    if sum(cond_vec)==length(cond_vec)
        
        %figure title
        % write cell's ID
        str = sprintf('Cell #%d, Date: %d, Bat %d',cell_num,day,bat) ;
        annotation('textbox',[.06 .88 .3 .1],'string',str,'EdgeColor','none','fontsize',15)
        %         for ii=1:length(k_sizes)
        %             for th_i=1:2
        %                 th=th_vec(th_i);
        for solo_thresh=1:2
            if solo_thresh==1
                solo_data=solo_25;
                str='Solo threshold = 25';
            else
               solo_data=solo_40;
                str='Solo threshold = 40';
            end
            for ii_dir=1:2
                
             
                raster_spike_pos=solo_data.solo(ii_dir).spikes.x_pos(:);
                raster_spike_ts=solo_data.solo(ii_dir).spikes.ts_usec(:);
                raster_bsp_pos=solo_data.solo(ii_dir).bsp.x_pos(:);
                raster_bsp_ts=solo_data.solo(ii_dir).bsp.ts_usec(:);
                time_limits = [min(solo_data.solo(ii_dir).bsp.ts_usec(:)) max(solo_data.solo(ii_dir).bsp.ts_usec(:))];
                bins=solo_data.solo(1).x_pos_firing_rate{1, 2};
                PSTH=solo_data.solo(ii_dir).PSTH_for_field_detection;
                fields_height=solo_data.solo(ii_dir).field_height;
                fields_center=solo_data.solo(ii_dir).field_center;
                if length(solo_data.solo(ii_dir).fields)==0
                    fields_edges=[];
                    FE_field_pass_num=[];
                    num_flights_with_spikes=[];
                else
                    fields_edges=[reshape([solo_data.solo(ii_dir).fields.edges_prc],2,length([solo_data.solo(ii_dir).fields.edges_prc])/2)'];                    
                    FE_field_pass_num=[solo_data.solo(ii_dir).fields.FE_field_pass_num];
                    num_flights_with_spikes=[solo_data.solo(ii_dir).fields.num_flights_with_spikes];
                end
                plot_pos_1=[x_pos(ii_dir), y_pos(solo_thresh),raster_size];
                plot_pos_2=[x_pos(ii_dir), y_pos(solo_thresh)+raster_size(2),PSTH_size];
                title_str=str;
                plot_solo_rater(raster_spike_pos,raster_spike_ts,time_limits,bins,PSTH,fields_height,fields_center,plot_pos_1,plot_pos_2,title_str,fields_edges,FE_field_pass_num,num_flights_with_spikes,raster_bsp_pos,raster_bsp_ts)
            end
            %             end
            %         end
        end
        % save figure
        figure_name = [figure_output_folder,'\solo_threshold_bat',num2str(bat),'_day_',num2str(day),'_cell_', num2str(cell_num),'.jpg'];
        saveas(gcf,figure_name)
        clf
    end
end
%%
function plot_solo_rater(raster_pos,raster_ts,time_limits,bins,PSTH,fields_height,fields_center,plot_pos_1,plot_pos_2,title_str,fields_edges,FE_field_pass_num,num_flights_with_spikes,raster_bsp_pos,raster_bsp_ts)
tunnel_limits=[0 135];
% raster of solo spikes
axes('position',plot_pos_1)
hold on
plot(raster_bsp_pos,raster_bsp_ts,'.','color',[.8 .8 .8],'markersize',2)
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

max_x = 1.5*ceil(max(PSTH));
for field_i=1:length(FE_field_pass_num)
    plot(fields_edges(field_i,:),(fields_height(field_i)/2)*ones(1,2),'color',[.5 .5 .5])
    text(mean(fields_edges(field_i)),(fields_height(field_i))*1.5,sprintf('%d/%d',num_flights_with_spikes(field_i),FE_field_pass_num(field_i)),'HorizontalAlignment','center')
end
max_x=max_x+1;
if max_x~=0 & ~isnan(max_x)
    set(gca,'xlim',tunnel_limits,'XTick',[],'ylim',[0 max_x],'ytick',max_x);
end
str = sprintf('Hz');
ylabel(str)
box off
title(title_str)
end

