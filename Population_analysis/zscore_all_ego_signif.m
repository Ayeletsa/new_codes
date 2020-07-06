
close all

clear
sort_by_max=1;
%% data:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_info=dir(dir_data);
analysis_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\population\';


co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct\';
figure_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\population\';

inclusion_dir='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\inclusion_cells_struct\';
inclusion_dir_info=dir(inclusion_dir);
inclusion_names={inclusion_dir_info.name};
inclusion_names([inclusion_dir_info.isdir])=[];
cell_ind=regexp(inclusion_names{2},'cell_');

cell_nums_inclusion=(cellfun(@(c) c(cell_ind+5:cell_ind+7),inclusion_names,'UniformOutput',false));

%% params for valid cells during CO
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
pop_param_file_name=fullfile(param_folder,'co_population_params.mat');
load(pop_param_file_name)
per_field_param_file_name=fullfile(param_folder,'per_field_params.mat');
load(per_field_param_file_name)
co_param_file_name=fullfile(param_folder,'co_params.mat');
load(co_param_file_name)

%% initialize:
per_field_cell_dir_count=0;
fr_count_all=0;
cell_dir_per_field_count=0;

all_per_field_valid_tuning={};
all_pos_width=[];
all_neg_width=[];
all_compound_width=[];
all_compound_pos_width=[];
all_compound_neg_width=[];
all_pos_rise_time=[];
all_neg_rise_time=[];
all_compound_pos_rise_time=[];
all_compound_neg_rise_time=[];

%% TO DO:
%check why there are large fields
%hist of per field sizes
% set minimum spikes for per field analysis
%%
file_names={dir_info.name};
file_names=file_names(find([dir_info.isdir]==0));

%% 1. run over cells:
for cell_i=1:length(file_names)
    % load data:
    load(fullfile(dir_data,file_names{cell_i}))
    cell_num=cell_co_solo_initial_analysis.exp_data.cell_num;
    %load cell's inclusion:
    Match=cellfun(@(a) find(contains(a,num2str(cell_num))) , cell_nums_inclusion, 'UniformOutput', 0);
    r=find(cellfun(@(c) ~isempty(c),Match));
    load(fullfile(inclusion_dir,inclusion_names{r}))
    
    
    for dir_i=1:2
        
       % basic cell data:
        bat=cell_co_solo_initial_analysis.exp_data.bat;
        day=cell_co_solo_initial_analysis.exp_data.day;
        cell_num=cell_co_solo_initial_analysis.exp_data.cell_num;
        %ego tuning of the cell:
        all_ego_tuning(cell_i,dir_i,:)=cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:);
        all_zscore_ego_tuning(cell_i,dir_i,:)=zscore(cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:));

        %% ego signif cell
 
%         shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
%         shuffle_file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);

        ego_signif_cells(cell_i,dir_i)=inclusion(dir_i).ego_cell;
        pyramidal_cells(cell_i,dir_i)=inclusion(dir_i).pyr;
        place_cells(cell_i,dir_i)=inclusion(dir_i).place_cell;
        valid_cells(cell_i,dir_i)=inclusion(dir_i).valid_cell;
    end

end
   
%% plot
figure('units','normalized','outerposition',[0 0 1 1])
mat_size=[0.1 0.1];
tuning_size=[0.07 0.06];

hor_dist=0.15;
x_pos(1)=0.03;
x_pos(2)=x_pos(1)+hor_dist;
x_pos(3)=x_pos(2)+hor_dist;
x_pos(4)=x_pos(3)+hor_dist;
x_pos(5)=x_pos(4)+hor_dist;
x_pos(6)=x_pos(5)+hor_dist;

ver_dist=0.2;
y_pos_mat(1)=0.7;
%% plot
ylimits=[-1 1];

% tuning all ego cells
ix=1;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=[squeeze(all_zscore_ego_tuning(ego_signif_cells(:,1),1,:));squeeze(all_zscore_ego_tuning(ego_signif_cells(:,2),2,:))];
title_str='zscore tuning of ego signif';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_X_bins_vector_of_centers,title_str,sort_by_max,ylimits)

% pyramidal ego
pyramidal_ego=ego_signif_cells&pyramidal_cells;
ix=2;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=[squeeze(all_zscore_ego_tuning(pyramidal_ego(:,1),1,:));squeeze(all_zscore_ego_tuning(pyramidal_ego(:,2),2,:))];
title_str='zscore tuning of pyramidal ego signif';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_X_bins_vector_of_centers,title_str,sort_by_max,ylimits)

% interneurons ego
non_pyramidal_ego=ego_signif_cells&~pyramidal_cells;
ix=3;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=[squeeze(all_zscore_ego_tuning(non_pyramidal_ego(:,1),1,:));squeeze(all_zscore_ego_tuning(non_pyramidal_ego(:,2),2,:))];
title_str='zscore tuning of interneurons ego signif';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_X_bins_vector_of_centers,title_str,sort_by_max,ylimits)

% all interneurons 
non_pyramidal=valid_cells&~pyramidal_cells;
ix=4;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=[squeeze(all_zscore_ego_tuning(non_pyramidal(:,1),1,:));squeeze(all_zscore_ego_tuning(non_pyramidal(:,2),2,:))];
data(nansum(data,2)==0,:)=[];
title_str='zscore tuning of interneurons';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_X_bins_vector_of_centers,title_str,sort_by_max,ylimits)

% all cells 
ix=5;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=[squeeze(all_zscore_ego_tuning(:,1,:));squeeze(all_zscore_ego_tuning(:,2,:))];
data(nansum(data,2)==0,:)=[];
title_str='zscore tuning ALL cells';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_X_bins_vector_of_centers,title_str,sort_by_max,ylimits)

% place cells 
ix=6;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=[squeeze(all_zscore_ego_tuning(place_cells(:,1),1,:));squeeze(all_zscore_ego_tuning(place_cells(:,2),2,:))];
data(nansum(data,2)==0,:)=[];
title_str='zscore tuning of place cells';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_X_bins_vector_of_centers,title_str,sort_by_max,ylimits)

%% save
fig_file_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\population\zscore_all_ego_signif',num2str(sort_by_max),'.jpg'];
saveas(gcf,fig_file_name)

%%
function plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)
axes('position',pos_ax_tuning)
plot(dis_per_field_bin_vec_of_center,nanmean(data))
title(title_str)
ylim(ylimits)
xlim([-40 40])
set(gca,'xtick',[])


normalize_to_other_max=[];
%[xlimits, ylimits] = fn_plot_2D_field (data, dis_per_field_X_bins_vector,dis_per_field_bin_vec_of_center, cells_bin,cells_bin_vector_of_centers,normalize_to_other_max);
axes('position',pos_ax_mat)
if sort_by_max==1
[~,id]=max(data,[],2);
else
  [~,id]=min(data,[],2);
end
  
[~,idd]=sort(id);
data=data(idd,:);
imagesc(data)
set(gca,'tickdir','out')

xlabel('Inter-bat distance (m)')
set(gca,'xtick',0:10:length(dis_per_field_bin_vec_of_center), 'xticklabel',[-40 -20 0 20 40])
%ylabel('Fields')
colormap(jet)
colorbar('location','eastoutside')


end
