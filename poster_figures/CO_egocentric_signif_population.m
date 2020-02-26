close all
clear
clc
all=1;
%fig_folder_name='D:\Ayelet\2bat_proj\Analysis\co_analysis\dis_bin_3_allo_bin_3';
%out_dir_type='D:\Ayelet\2bat_proj\Analysis\co_analysis\ego_signif_cells';

out_pop_fig_dir='D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures';

unstabel_cells=[27, 36, 228,263];
%inhibit_cells=[215,2;218,2;233,2;251,1;254,1;257,1;267,2;270,2;287,2;];
inhibit_cells=[8, 17, 21, 31, 32, 53, 57, 74,215,218,233,251,254,257,267,270,287];
%% structs' parameters

behavior_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
files = dir(behavior_structs_folder);
behavior_struct_names = {files.name};
co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct';
files = dir(co_shuffle_structs_folder);
co_shuffle_struct_names = {files.name};

%% params for valid cells during CO

min_spikes = 30;
min_ego_inf = 0.1;
min_even_odd = 0.2;
alpha = .05;
interneuron_firing_rate = 10;
max_for_pyramidal=5;
%% figure prop


figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [24 21]; % ~A4 size
set(gcf,'DefaultAxesFontSize',9);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');
line_width=3;
ego_color=[1  0 1];
fsize = 12;

% tuning cuvre mat over all cells
ax_plot{1}=axes('position',[0.08 0.55 0.2 0.2]); %sorted by peak
ax_plot{2}=axes('position',[0.08 0.75 0.2 0.1]);%average over all cells


% z score:
ax_plot{3}=axes('Position',[0.43 0.55 0.2 0.2]); %mat
ax_plot{4}=axes('Position',[0.43 0.75 0.2 0.1]); %average

%%
ii_fr=1;
for ii_cell = 3:length(behavior_struct_names)
    
    %% load data
    
    struct_name =behavior_struct_names{ii_cell};
    file_name = fullfile(behavior_structs_folder,struct_name);
    load(file_name);
    if ii_cell<=102
    behavior_struct=cell_co_solo_initial_analysis;
    else
            load(file_name);

    end
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    cell_num=behavior_struct.exp_data.cell_num;
    if ismember(cell_num,unstabel_cells)
        continue
    else
        
        shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
        file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);
        if exist(file_name)
            co_shuffle_struct=load(file_name);
            co_shuffle_struct=co_shuffle_struct.shuffling_struct;
        end
        % loop for each direction
        for ii_dir = 1:2
            
            n_spikes = behavior_struct.co(ii_dir).info.n_spikes;
            even_odd_coherence = co_shuffle_struct(ii_dir).odd_even_coherence.corr;
            ego_inf = co_shuffle_struct(ii_dir).shuffled_data.params.information_per_spike_ego.values(1);
            ego_inf_p = co_shuffle_struct(ii_dir).shuffled_data.params.information_per_spike_ego.p_value;
            cv_p=co_shuffle_struct(ii_dir).shuffled_data.params.cv.p_value;
            
            % conditions :
            a = n_spikes > min_spikes;
            b = even_odd_coherence > min_even_odd;
            c = ego_inf > min_ego_inf;
            d = ego_inf_p < alpha;
            e=behavior_struct.exp_data.mean_fr<max_for_pyramidal;
            
            if min([a,b,c,d,e]) == 1
                % save normalized firing rate
                fr = co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate(1,:);
                zscore_tuning=co_shuffle_struct(ii_dir).shuffled_data.zscore_of_tuning_curve_data ;
                %                 norm_fr = fr/max(fr);
                range_fr = fr - min(fr);
                norm_fr = range_fr/max(range_fr);
                firing_rate_mat(ii_fr,:) = norm_fr;
                zscore_mat(ii_fr,:)=zscore_tuning;
                % save the cell's number and direction
                cells_num_dir(ii_fr,:) = [cell_num,ii_dir];
                
                ii_fr = ii_fr+1;
            end
        end
    end
end

%% Plot
bins=behavior_struct.co(1).firing_rate.dis_m(2,:) ;
inhibit_cells_ind=find(sum(cells_num_dir(:,1)==inhibit_cells,2));

%% 1. all cells tuning

%a. mat of all tunings for all cells normalized:
axes(ax_plot{1})
exs_mat=firing_rate_mat;
if all==0
    exs_mat(inhibit_cells_ind,:)=[];

end
    
mat=exs_mat;
[~, max_mat_ind]=max(mat,[],2);
[~, sort_max_mat_ind]=sort(max_mat_ind);
sorted_mat=mat(sort_max_mat_ind,:);
c_max=[];
plot_2d_mat(bins,1:size(sorted_mat,1),sorted_mat,c_max)
xlabel('Inter-bat distance (m)','fontsize',fsize,'fontweight','bold')
set(gca,'ytick',0:10:size(sorted_mat,1),'TickDir','out','TickLength',[0.02 0])
ylabel('Cells','fontsize',fsize,'fontweight','bold')
set(gca,'xtick',-40:20:40,'TickDir','out','TickLength',[0.02 0])

%b. average over all cells:
axes(ax_plot{2})
plot(bins, nanmean(sorted_mat),'lineWidth', line_width,'color',ego_color)
xlim([min(bins) max(bins)])
ylim([min(nanmean(sorted_mat)) max(nanmean(sorted_mat))*1.05]);
set(gca,'xtick',[],'TickDir','out','TickLength',[0.02 0])
ylabel(sprintf('Average\nover cells'),'fontsize',fsize,'fontweight','bold')
title(sprintf('Normalized firing rates'),'fontsize',fsize,'fontweight','bold')
box off
%c. mat of z scores all cells:
axes(ax_plot{3})
 exs_mat=zscore_mat;
 if all==0

exs_mat(inhibit_cells_ind,:)=[];
 end
 mat=exs_mat;mat=mat-min(mat(:));
[~, max_mat_ind]=max(mat,[],2);
[~, sort_max_mat_ind]=sort(max_mat_ind);
sorted_mat=mat(sort_max_mat_ind,:);
c_max=[];
plot_2d_mat(bins,1:size(sorted_mat,1),sorted_mat,c_max)
xlabel('Inter-bat distance (m)','fontsize',fsize,'fontweight','bold')
ylabel('Cells','fontsize',fsize,'fontweight','bold')
set(gca,'ytick',0:10:size(sorted_mat,1),'TickDir','out','TickLength',[0.02 0])
set(gca,'xtick',-40:20:40,'TickDir','out','TickLength',[0.02 0])

%d. average zscore:
axes(ax_plot{4})
plot(bins, nanmean(zscore_mat),'lineWidth', line_width,'color',ego_color)
xlim([min(bins) max(bins)])
ylim([min(nanmean(zscore_mat)) max(nanmean(zscore_mat))*1.05]);
set(gca,'xtick',[],'TickDir','out','TickLength',[0.02 0])
ylabel(sprintf('Average\nover cells'),'fontsize',fsize,'fontweight','bold')
title(sprintf('Z-scored firing rates'),'fontsize',fsize,'fontweight','bold')
box off

%% save
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
if all==0
    fig_name=fullfile(dir_name,'population_egocentric_cells_enhanced');

else
fig_name=fullfile(dir_name,'population_egocentric_cells_all');
end
print(gcf, fig_name, '-dpdf',  '-painters');
