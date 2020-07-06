%% Population per field z score
clear
%general params:
sort_by_max=0;

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
per_field_param_file_name=fullfile(param_folder,'per_field_params.mat');
load(per_field_param_file_name)

%% initialize:
per_field_cell_dir_count=0;
fr_count_all=0;
cell_dir_per_field_count=0;

all_per_field_valid_tuning=[];
all_per_field_valid_tuning_zscore=[];
all_pos_width=[];
all_neg_width=[];
all_compound_width=[];
all_compound_pos_width=[];
all_compound_neg_width=[];
all_pos_rise_time=[];
all_neg_rise_time=[];
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
       
        %% ego signif cell
 
%         shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
%         shuffle_file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);

        ego_signif_cells(cell_i,dir_i)=inclusion(dir_i).ego_cell;
        
        %% per field data
        %1. Test if run on per field:
       
%         a=sum(~isnan(cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(dir_i).spikes.ts_usec(:)))>=min_n_spike;
%         b=cell_co_solo_initial_analysis.solo(dir_i).SI>SI_threshold;
%         c=~isempty(cell_co_solo_initial_analysis.solo(dir_i).fields);
%         d=cell_co_solo_initial_analysis.exp_data.mean_fr<max_for_pyramidal;
%         per_field_cell_cond=a & b & c & d;
        if inclusion(dir_i).place_cell 
          cell_dir_per_field_count=cell_dir_per_field_count+1;

            per_field_data=cell_co_solo_initial_analysis.co(dir_i).per_field_href;
            solo_field_data=cell_co_solo_initial_analysis.solo(dir_i).fields;
            
          

          for field_i=1:length(solo_field_data)
              % run only on per fields that obey conditions:
              r=per_field_data(field_i).tuning_dis_x_fr_per_field;
              [ind_length,~,~]=find_length_of_consecutive_ind(find(~isnan(r)),length(r));
              a=max(ind_length)>=min_r_length_per_field*length(r); %enough data that is not nan in tuning curve
              b=per_field_data(field_i).number_of_spikes_per_field>=min_n_spike_per_field;% enough spikes
              per_field_cond=a & b;
               if per_field_cond
                fr_count_all=fr_count_all+1;
                %1. per field tuning:
                all_per_field_valid_tuning(fr_count_all,:)=nan*zeros(size(r));
                all_per_field_valid_tuning(fr_count_all,:)=r';
                non_nan_ind=~isnan(r);
                all_per_field_valid_tuning_zscore(fr_count_all,:)=nan*zeros(size(r));
                all_per_field_valid_tuning_zscore(fr_count_all,non_nan_ind)=zscore(r(non_nan_ind))';
                
                % solo_fe
%                 % iso dis:
%                 Isolation_dis(cell_dir_per_field_count)=cell_co_solo_initial_analysis.exp_data.Isolation_dis;
% 
%                 %2. get solo prop per field:
                 per_field_solo_height(fr_count_all)=solo_field_data(field_i).peak;
%                 per_field_solo_nrm_height(fr_count_all)=solo_field_data(field_i).peak./nanmean(cell_co_solo_initial_analysis.solo(dir_i).PSTH_for_field_detection);
%                 per_field_solo_size(fr_count_all)=solo_field_data(field_i).width_href;
%                 per_field_solo_SI(fr_count_all)=solo_field_data(field_i).field_SI;
%                  
%                 %3. get co prop per field:
%                 per_field_co_SI(fr_count_all)=per_field_data(field_i).SI;
%                 per_field_co_CV(fr_count_all)=per_field_data(field_i).cv;
%                 per_field_co_modulation_depth(fr_count_all)=per_field_data(field_i).modulation_depth;
%                 per_field_co_sparsity(fr_count_all)=per_field_data(field_i).sparsity;
%                 
                %4. get tuning width and rise time:
                if per_field_data(field_i).time_signif_field.signif_based_on_extreme_bins==1 
                    rise_and_width_data=get_per_field_tuning_width_and_rise_time(per_field_data(field_i).time_signif_field,min_dis_pos_neg);
                    id_per_field_non_signif(fr_count_all)=0;
                    id_per_field_pos_tuning(fr_count_all)=rise_and_width_data.pos_per_field_tuning;
                    id_per_field_neg_tuning(fr_count_all)=rise_and_width_data.neg_per_field_tuning;
                    
% 
%                     all_pos_width=[all_pos_width,rise_and_width_data.pos_width];
%                     all_neg_width=[all_neg_width,rise_and_width_data.neg_width];
%                     all_compound_width=[all_compound_width,rise_and_width_data.compound_width];
%                     all_compound_pos_width=[all_compound_pos_width,rise_and_width_data.pos_compound_width];
%                     all_compound_neg_width=[all_compound_neg_width,rise_and_width_data.neg_compound_width];
%                     all_pos_rise_time=[all_pos_rise_time,rise_and_width_data.pos_rise_time];
%                     all_neg_rise_time=[all_neg_rise_time,rise_and_width_data.neg_rise_time];
                    
              
                else
                    id_per_field_non_signif(fr_count_all)=1;
                    id_per_field_pos_tuning(fr_count_all)=0;
                    id_per_field_neg_tuning(fr_count_all)=0;
                  
                end
               end
          end
%           % compute cell per field corr:
%           if size(all_per_field_valid_tuning,1)>1
%               per_field_cell_dir_count=per_field_cell_dir_count+1;
%               corrs=tril(corr(all_per_field_valid_tuning',all_per_field_valid_tuning','rows','pairwise'),-1);
%               corr_per_field_per_cell(per_field_cell_dir_count)=nanmean(corrs(corrs~=0));
%           
%           end
        end
    end
%     %% compute corr between dirs:
%     corrs=tril(corr(squeeze(all_ego_tuning(cell_i,:,:))','rows','pairwise'),-1);
%     corr_ego_tuning_between_dirs(cell_i)=corrs(2,1);
end
   %% divde solo fr to 4 subdevisions:
   marings=prctile(per_field_solo_height,[25,50,75]);
   field_group{1}=find(per_field_solo_height<marings(1));
   field_group{2}=find(per_field_solo_height>marings(1) & per_field_solo_height<marings(2));
   field_group{3}=find(per_field_solo_height>marings(2) & per_field_solo_height<marings(3));
   field_group{4}=find(per_field_solo_height>marings(3));

%% plot
figure('units','normalized','outerposition',[0 0 1 1])
mat_size=[0.08 0.08];
tuning_size=[0.05 0.06];

hor_dist=0.095;
x_pos(1)=0.03;
x_pos(2)=x_pos(1)+hor_dist;
x_pos(3)=x_pos(2)+hor_dist;
x_pos(4)=x_pos(3)+hor_dist;
x_pos(5)=x_pos(4)+hor_dist;
x_pos(6)=x_pos(5)+hor_dist;
x_pos(7)=x_pos(6)+hor_dist;
x_pos(8)=x_pos(7)+hor_dist;
x_pos(9)=x_pos(8)+hor_dist;
x_pos(10)=x_pos(9)+hor_dist;

ver_dist=0.2;
y_pos_mat(1)=0.83;
y_pos_mat(2)=y_pos_mat(1)-ver_dist;
y_pos_mat(3)=y_pos_mat(2)-ver_dist;
y_pos_mat(4)=y_pos_mat(3)-ver_dist;
y_pos_mat(5)=y_pos_mat(4)-ver_dist;
%%
ylim_tuning=[0 15];
y_lim_zscore=[-1 1];
% tuning:
ylimits=ylim_tuning;
ix=1;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning;
title_str='tuning of all valid fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

ix=1;
iy=2;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(~id_per_field_non_signif,:);
title_str='tuning of signif fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=1;
iy=3;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(id_per_field_pos_tuning,:);
title_str='tuning of pos signif fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=1;
iy=4;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(id_per_field_neg_tuning,:);
title_str='tuning of neg signif fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


id_compound=intersect(find(id_per_field_neg_tuning),find(id_per_field_pos_tuning));
ix=1;
iy=5;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(id_compound,:);
title_str='tuning of compound fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


%zcore
ylimits=y_lim_zscore;
ix=2;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore;
title_str='zscore tuning of all valid fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

ix=2;
iy=2;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(~id_per_field_non_signif,:);
title_str='zscore tuning of signif fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

ix=2;
iy=3;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(id_per_field_pos_tuning,:);
title_str='zscore tuning of pos signif fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=2;
iy=4;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(id_per_field_neg_tuning,:);
title_str='zscore tuning of neg signif fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=2;
iy=5;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(id_compound,:);
title_str='zscore tuning of compound fields';
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

%% for sub groups of solo fields
for i=1:4
id_cells=field_group{i};
%tuning:
ylimits=ylim_tuning;

ix=2*i+1;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(id_cells,:);
title_str=sprintf('tuning all group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

ix=2*i+1;
iy=2;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(intersect(find(~id_per_field_non_signif),id_cells),:);
title_str=sprintf('tuning signif group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=2*i+1;
iy=3;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(intersect(find(id_per_field_pos_tuning),id_cells),:);
title_str=sprintf('tuning pos group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=2*i+1;
iy=4;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(intersect(find(id_per_field_neg_tuning),id_cells),:);
title_str=sprintf('tuning neg group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


id_compound=intersect(find(id_per_field_neg_tuning),find(id_per_field_pos_tuning));
ix=2*i+1;
iy=5;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning(intersect(find(id_compound),id_cells),:);
title_str=sprintf('tuning compound group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

%zscore
ylimits=y_lim_zscore;

ix=2*i+2;
iy=1;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(id_cells,:);
title_str=sprintf('zscore all group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

ix=2*i+2;
iy=2;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(intersect(find(~id_per_field_non_signif),id_cells),:);
title_str=sprintf('zscore signif group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=2*i+2;
iy=3;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(intersect(find(id_per_field_pos_tuning),id_cells),:);
title_str=sprintf('zscore pos group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)


ix=2*i+2;
iy=4;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
data=all_per_field_valid_tuning_zscore(intersect(find(id_per_field_neg_tuning),id_cells),:);
title_str=sprintf('zscore neg group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

ix=2*i+2;
iy=5;
pos_ax_mat=[x_pos(ix),y_pos_mat(iy),mat_size];
pos_ax_tuning=[x_pos(ix),y_pos_mat(iy)+mat_size(2),tuning_size];
id_compound=intersect(find(id_per_field_neg_tuning),find(id_per_field_pos_tuning));
data=all_per_field_valid_tuning_zscore(intersect(find(id_compound),id_cells),:);
title_str=sprintf('zscore compound group %d',i);
plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)

end

%% save
fig_file_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\population\per_field_zscore_color_bar_sort_by_',num2str(sort_by_max),'.jpg'];
saveas(gcf,fig_file_name)


%%
function plot2dmats(pos_ax_mat,pos_ax_tuning, data, dis_per_field_bin_vec_of_center,title_str,sort_by_max,ylimits)
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
xlabel('Inter-bat distance (m)')
set(gca,'xtick',0:10:length(dis_per_field_bin_vec_of_center), 'xticklabel',[-40 -20 20 40])
%ylabel('Fields')
colormap(jet)
colorbar('location','eastoutside')

axes('position',pos_ax_tuning)
plot(dis_per_field_bin_vec_of_center,nanmean(data))
title(title_str)
ylim(ylimits)
xlim([-40 40])
set(gca,'xtick',[])
end
