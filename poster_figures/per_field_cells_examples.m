
% plot imporatnat statistics for cross-overs analysis
close all
clear
clc
cells_num_dir=[294 2;257 2;261 2;277 1; 277 2;282 1; 296 2;267 2];
relevant_field{1}=[1,2,3];
relevant_field{2}=[1,2];
relevant_field{3}=[1,2];
relevant_field{4}=[1,2,3];
relevant_field{5}=[1,2];
relevant_field{6}=[1,2];
relevant_field{7}=[1];
relevant_field{8}=[1,2,3];






%% structs' parameters

CO_solo_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
files = dir(CO_solo_structs_folder);
behavior_struct_names = {files.name};
co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct\';
files = dir(co_shuffle_structs_folder);
co_shuffle_struct_names = {files.name};


%% figure parameters

spike_colors = {[1  0 1],[1  0 1],[.95 .95 .95]};
solo_colors = {[0 0 0],[0 0 0]};
FR_colors = {[1 1 1],[1 1 1]};%
coherence_colors = {[1 .5 0],[.8 0 .6]};
tunnel_limits = [0 122];
sig_bkg_color = [.8 .95 1];
lwidth = 2;
alpha = .05;


%% open figure

figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [24 21]; % ~A4 size
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');

fsize = 10;

horizontal_dis=0.45;
x_main(1)=0.08;
x_main(2)=x_main(1)+horizontal_dis;
x_main(3:2:8)=x_main(1);
x_main(4:2:8)=x_main(2);


x_ego=x_main;
x_allo=x_main+0.13;
x_pos_per_field=x_allo+0.06;

vert_dis=0.24;
y_main(1:2)=0.8;
y_main(3:4)=y_main(1)-vert_dis;
y_main(5:6)=y_main(3)-vert_dis;
y_main(7:8)=y_main(5)-vert_dis;


y_ego=y_main+0.12;
y_allo=y_main;
y_pos_per_field=y_main;

for ii_cell = 1:length(cells_num_dir)
              
    ego_tun_ax{ii_cell} = axes('units','normalized','Position',[x_ego(ii_cell) y_ego(ii_cell) 0.12 0.05]);
    
    main_ax{ii_cell} = axes('units','normalized','Position',[x_main(ii_cell) y_main(ii_cell) 0.12 0.12]); %
    
    allo_ax{ii_cell}=axes('units','normalized','Position',[x_allo(ii_cell) y_allo(ii_cell) 0.03 0.12]);
end



%% plot for each cell

for ii_cell = 1:length(cells_num_dir)
    fields_to_plot=relevant_field{ii_cell};
    num_fileds=length(fields_to_plot);
    bin_dis_i=1;
    bin_allo_i=1;
    %% load data
    cell_file_ind=find(~cellfun(@isempty,regexp(behavior_struct_names,sprintf('cell_%d.mat',cells_num_dir(ii_cell,1)))));
    ii_dir=cells_num_dir(ii_cell,2);
    
    struct_name =behavior_struct_names{cell_file_ind};
    file_name = fullfile(CO_solo_structs_folder,struct_name);
    load(file_name);
    behavior_struct=cell_co_solo_initial_analysis;
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    cell_num=behavior_struct.exp_data.cell_num;
    
    shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
    file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);
    if exist(file_name)
        co_shuffle_struct=load(file_name);
        co_shuffle_struct=co_shuffle_struct.shuffling_struct;
        
        dis_before_after_co = behavior_struct.co(1).params.dis_before_after_co;
        %%
        axes(main_ax{ii_cell})
        hold on
        %         line([0 0],tunnel_limits,'Color',[.7 .7 .7],'LineStyle','--')
        %         plot(behavior_struct.co(ii_dir).bsp.dis_m(:),behavior_struct.co(ii_dir).bsp.x_pos(:),'.','color',spike_colors{3})
        %         plot(behavior_struct.co(ii_dir).spikes.dis_m(:),behavior_struct.co(ii_dir).spikes.x_pos(:),'.','color',spike_colors{ii_dir},'markersize',8)
        normalize_to_other_max=[];
        field_density_mat_X_Y=behavior_struct.co(ii_dir).firing_rate.field_density_smoothed_XY_with_NaN{bin_dis_i,bin_allo_i};
        X_bins_vector=behavior_struct.co(ii_dir).firing_rate.dis_X_bins_vector{bin_dis_i,bin_allo_i};
        X_bins_vector_of_centers=behavior_struct.co(ii_dir).firing_rate.dis_X_bins_vector_of_centers{bin_dis_i,bin_allo_i};
        Y_bins_vector=behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector{bin_dis_i,bin_allo_i};
        Y_bins_vector_of_centers=behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers{bin_dis_i,bin_allo_i};
        
        allo_bin_size=behavior_struct.co(ii_dir).firing_rate.allo_bin_size{bin_dis_i,bin_allo_i};
        dis_bin_size=behavior_struct.co(ii_dir).firing_rate.dis_bin_size{bin_dis_i,bin_allo_i};
        
        [xlimits, ylimits] = fn_plot_2D_field (field_density_mat_X_Y, X_bins_vector, X_bins_vector_of_centers,Y_bins_vector, Y_bins_vector_of_centers,normalize_to_other_max);
        set(gca,'xlim',xlimits,...
            'ylim',tunnel_limits,...
            'ytick',tunnel_limits)
        if ii_cell==7
        ylabel('Position in tunnel (m)','fontsize',fsize)
        xlabel('Inter-bat distance (m)','fontsize',fsize)
        end
        
        
        max_fr=max(field_density_mat_X_Y(:));
        text(xlimits(1)+xlimits(1),0.93*ylimits(2),sprintf('%.1f Hz',max_fr))
        box off
        
        %% egocentric firing rate
        axes(ego_tun_ax{ii_cell})
        hold on
        % plot(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers,co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate(2:end,:),'color',spike_colors{3},'LineWidth',lwidth);
        plot(behavior_struct.co(ii_dir).firing_rate.dis_X_bins_vector_of_centers{bin_dis_i,bin_allo_i}  ,behavior_struct.co(ii_dir).firing_rate.dis_x_pos_fr_for_2D{bin_dis_i,bin_allo_i},'color',spike_colors{ii_dir},'LineWidth',lwidth);
        max_y = round(max(max(behavior_struct.co(ii_dir).firing_rate.dis_x_pos_fr_for_2D{bin_dis_i,bin_allo_i}))*1.1*10) / 10;
        if max_y < 1
            max_y = 1;
        end
        bins_diff = mean(diff(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers));
        x_limits = [co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(1)-bins_diff/2 co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(end)+bins_diff/2];
        set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
        if ii_cell==7
        str = sprintf('Hz');
        ylabel(str,'fontsize',fsize)
        end
        box off
        title(sprintf('Cell %d',cell_num),'fontsize',fsize)
        
        %% solo firing rate - regular bins
        axes(allo_ax{ii_cell})
        hold on
        %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
        plot(behavior_struct.solo(ii_dir).PSTH_for_field_detection,behavior_struct.solo(1).x_pos_firing_rate{1, 2}   ,'color',solo_colors{ii_dir},'LineWidth',2);
   
        plot(behavior_struct.solo(ii_dir).field_height(fields_to_plot),behavior_struct.solo(ii_dir).field_center(fields_to_plot),'*r')
        if ~isempty(fields_to_plot)
            for fr_i= fields_to_plot
                plot([-0.5 -0.5],[behavior_struct.solo(ii_dir).field_edges(1,fr_i),behavior_struct.solo(ii_dir).field_edges(2,fr_i)],'color',[.5 .5 .5],'linewidth',2)
            end
        end
        %plot(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D{bin_dis_i,bin_allo_i},behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers{bin_dis_i,bin_allo_i} ,'color',spike_colors{ii_dir},'LineWidth',lwidth);
        max_x = 1.1*ceil(max(behavior_struct.solo(ii_dir).PSTH_for_field_detection));
        if max_x~=0 & ~isnan(max_x)
            set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[-1 max_x],'Xtick',max_x,'color',FR_colors{2});
        end
        plot([0 0], tunnel_limits,'k')
        if ii_cell==7
        str = sprintf('Hz');
        xlabel(str,'fontsize',fsize)
        end
        axis off
        
        % per field firing rate:
        if num_fileds>0
           % sorted_hights_ind_sorted=sort(sorted_hights_ind);
            dis_y_pos=0.04;
            y_pos(1)=y_pos_per_field(ii_cell);
            y_pos(2)=y_pos(1)+dis_y_pos;
            y_pos(3)=y_pos(2)+dis_y_pos;
            y_pos(4)=y_pos(3)+dis_y_pos;
            count=0;
            for fr_i=sort(fields_to_plot)
                count=count+1;
                
                axes('units','normalized','Position',[x_pos_per_field(ii_cell) y_pos(count) 0.1 0.03]);
                if sum(~isnan(behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1, fr_i}))~=0
                    bins_center=behavior_struct.co(ii_dir).firing_rate.dis_X_bins_vector_of_centers{bin_dis_i,bin_allo_i};
                    shuf_data=behavior_struct.co(ii_dir).firing_rate.signif_field{1, fr_i}.shuffled_data;
                    r=behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1, fr_i};
                    plot(bins_center,shuf_data,'color',spike_colors{3},'LineWidth',lwidth); hold on;
                    plot(bins_center ,r,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    
                    if behavior_struct.co(ii_dir).firing_rate.signif_field{1, fr_i}.signif_based_on_extreme_bins==1
                        neg_signif=behavior_struct.co(ii_dir).firing_rate.signif_field{1, fr_i}.neg_signif;
                        plot(bins_center(neg_signif),r(neg_signif)+1,'r*')
                        pos_signif=behavior_struct.co(ii_dir).firing_rate.signif_field{1, fr_i}.pos_signif;
                        plot(bins_center(pos_signif),r(pos_signif)+1,'g*')
                    end
                    max_y = round(max(max([r;shuf_data]))*1.1*10) / 10 +1;
                    if max_y < 1
                        max_y = 1;
                    end
                    bins_diff = mean(diff(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers));
                    x_limits = [co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(1)-bins_diff/2 co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(end)+bins_diff/2];
                    if count==1
                        set(gca,'xlim',x_limits,'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                    else
                        set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                        
                    end
                    r=behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1, fr_i};
%                     sparsity=(nanmean(r).^2)/nanmean(r.^2);
%                     std_r=nanstd(r);
%                     mean_r=nanmean(r);
%                     cv=std_r/mean_r;
%                     if behavior_struct.co(ii_dir).firing_rate.signif_field{1, fr_i}.signif_based_on_extreme_bins==1
%                         text(x_limits(2),max_y*0.5,sprintf('SI=%.2f sig bins \n #spikes=%d',behavior_struct.co(ii_dir).firing_rate.information_per_spike_per_field{fr_i},behavior_struct.co(ii_dir).firing_rate.number_of_spikes_per_field{fr_i}),'color',[1 0 0])
%                     else
%                         text(x_limits(2),max_y*0.5,sprintf('SI=%.2f \n #spikes=%d',behavior_struct.co(ii_dir).firing_rate.information_per_spike_per_field{fr_i},behavior_struct.co(ii_dir).firing_rate.number_of_spikes_per_field{fr_i}))
%                         
%                     end
%                     if behavior_struct.co(ii_dir).firing_rate.signif_field{1, fr_i}.signif_based_on_CV==1
%                         
%                         text(x_limits(2),0,sprintf('cv=%.2f',cv),'color',[1 0 0])
%                     else
%                         text(x_limits(2),0,sprintf('cv=%.2f',cv))
%                     end
%                     
%                     % str = sprintf('Hz');
                    %ylabel(str,'fontsize',fsize)
                    box off
                end
            end
        end
        
        
    end
end
%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,'per_field_cells');
print(gcf, fig_name, '-dpdf',  '-painters');

clf


