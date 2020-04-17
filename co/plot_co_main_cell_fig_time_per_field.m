function plot_co_main_cell_fig_time_per_field(dir_param_file_name,population_param_file_name,co_param_file_name)
% plot imporatnat statistics for cross-overs analysis
load(dir_param_file_name)
load(population_param_file_name)
load(co_param_file_name)
%% structs' parameters

files = dir(cell_co_solo_initial_analysis_struct_folder);
behavior_struct_names = {files.name};
files = dir(co_shuffle_folder_name);
co_shuffle_struct_names = {files.name};


%% figure parameters

spike_colors = {[1  0 1],[1  0 1],[.95 .95 .95]};
solo_colors = {[0 0 0],[0 0 0]};
FR_colors = {[1 1 1],[1 1 1]};%
coherence_colors = {[1 .5 0],[.8 0 .6]};
tunnel_limits = [0 135];
sig_bkg_color = [.8 .95 1];
lwidth = 3;
alpha = .05;


%% open figure

figure('units','normalized','outerposition',[0 0 1 1],'color',[1 1 1])
set(gcf,'DefaultAxesFontSize',12);
fsize = 14;


%% plot for each cell
for ii_cell = 162:length(behavior_struct_names)
    ii_cell
    signif=0;
    %% load data
    
    struct_name =behavior_struct_names{ii_cell};
    file_name = fullfile(cell_co_solo_initial_analysis_struct_folder,struct_name);
    load(file_name);
    behavior_struct=cell_co_solo_initial_analysis;
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    cell_num=behavior_struct.exp_data.cell_num;
    if ~isnan(behavior_struct.co(1).inter_field_firing_rate.inter_field_edge)  | ~isnan(behavior_struct.co(2).inter_field_firing_rate.inter_field_edge)
        shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
        file_name = fullfile(co_shuffle_folder_name,shuffle_struct_name);
        if exist(file_name)
            co_shuffle_struct=load(file_name);
            co_shuffle_struct=co_shuffle_struct.shuffling_struct;
            
            
            %% write cell's ID and stability
            % different bin sixe for 2D
            for bin_dis_i=1:size(behavior_struct.co(1).firing_rate.dis_bin_size,1)
                %  for bin_allo_i=1:size(behavior_struct.co(1).firing_rate.allo_bin_size,2)
                bin_allo_i=bin_dis_i;
                dis_before_after_co = behavior_struct.co(1).params.dis_before_after_co;
                
                % write cell's ID
                str = sprintf('Cell #%d, Date: %d, Bat %d',cell_num,day,bat) ;
                annotation('textbox',[.06 .88 .3 .1],'string',str,'EdgeColor','none','fontsize',15)
                
                % cell stability
                n_bahviors = numel(behavior_struct.exp_data.stability{1, 1});
                sleep_ind = [1 n_bahviors];
                behaving_ind = 2:n_bahviors-1;
                axes('units','normalized','Position',[.07 .87 .1 .05]);
                min_t = min(cell2mat(behavior_struct.exp_data.stability{1, 1})) - mean(diff(behavior_struct.exp_data.stability{1, 1}{1,1}));
                bar(cell2mat(behavior_struct.exp_data.stability{1, 1}(1,behaving_ind))-min_t,cell2mat(behavior_struct.exp_data.stability{1, 2}(1,behaving_ind)),1,'facecolor',[.7 .7 .7])
                hold on
                bar(cell2mat(behavior_struct.exp_data.stability{1, 1}(1,sleep_ind))-min_t,cell2mat(behavior_struct.exp_data.stability{1, 2}(1,sleep_ind)),1,'facecolor',[0 .5 1])
                max_fr_in_bin = max(cell2mat(behavior_struct.exp_data.stability{1, 2}));
                max_t = max(cell2mat(behavior_struct.exp_data.stability{1, 1})) - min_t;
                max_t_minutes = max_t/(1e6*60);
                if max_fr_in_bin>0
                    set(gca,'xlim',[0 max_t*1.02],'xtick',[0 max_t*1.02],'xticklabel',[0 round(max_t_minutes*1.02)],...
                        'ylim',[0 max_fr_in_bin*1.1],'ytick',[0 max_fr_in_bin*1.1],'yticklabel',[0 round(max_fr_in_bin*1.1)])%,'TickLabelRotation ',45)
                else
                    set(gca,'xlim',[0 max_t*1.02],'xtick',[0 max_t*1.02],'xticklabel',[0 round(max_t_minutes*1.02)],...
                        'ylim',[0 1*1.1],'ytick',[0 1*1.1],'yticklabel',[0 round(1*1.1)])%,'TickLabelRotation ',45)
                end
                xlabel('Time (min)','fontsize',fsize*.7)
                ylabel('Hz')
                box off
                
                % ratio between pre and post sleeps
                
                axes('units','normalized','Position',[.19 .87 .02 .05]);
                stability_index = behavior_struct.exp_data.stability{1, 3}{1, 3};
                if ~isnan(stability_index)
                    bar(1,stability_index,.2,'facecolor',[0 .5 1])
                    set(gca,'xtick',1,'xticklabel','Pre/Post','ylim',[0 stability_index*1.05],'ytick',stability_index,'yticklabel',round(stability_index*10)/10)
                    box off
                end
                % Mean firing rate along the entire session
                axes('units','normalized','Position',[.23 .87 .02 .05]);
                bar(1,behavior_struct.exp_data.mean_fr  ,.2,'facecolor',[1 .5 0])
                set(gca,'xtick',1,'xticklabel','mean FR','ylim',[0 behavior_struct.exp_data.mean_fr  *1.05],'ytick',behavior_struct.exp_data.mean_fr  ,'yticklabel',round(behavior_struct.exp_data.mean_fr  *10)/10)
                box off
                
                %%     plot data for every direction
                
                for ii_dir = 1:2
                    
                    dir_adj = (ii_dir-1)*.5;
                    % main figure
                    co_ax{1} = axes('units','normalized','Position',[0.05+dir_adj 0.58 0.12 0.17]); % main fig
                    co_ax{2} = axes('units','normalized','Position',[0.05+dir_adj 0.75 0.12 0.06]); % ego fr
                    co_ax{3} = axes('units','normalized','Position',[0.17+dir_adj 0.58 0.03 0.17]); % allo fr
                    % allocentric fr comparison
                    co_ax{4} = axes('units','normalized','Position',[0.22+dir_adj 0.58 0.03 0.17]);
                    % solo raster plot
                    co_ax{5} = axes('units','normalized','Position',[0.27+dir_adj 0.58 0.03 0.17]); % allo fr comparison
                    co_ax{20} = axes('units','normalized','Position',[0.3+dir_adj 0.58 0.03 0.17]); % allo fr comparison
                    
                    % 2D representation
                    %                 co_ax{6} = axes('units','normalized','Position',[0.37+dir_adj 0.58 0.035 0.17]); % bsp only
                    %                 co_ax{7} = axes('units','normalized','Position',[0.40+dir_adj 0.58 0.035 0.17]); % during solo
                    %                 co_ax{8} = axes('units','normalized','Position',[0.43+dir_adj 0.58 0.035 0.17]); % during co
                    %                 co_ax{9} = axes('units','normalized','Position',[0.46+dir_adj 0.58 0.035 0.17]); % during solo+co
                    %                 % info
                    % co_ax{10} = axes('units','normalized','Position',[0.05+dir_adj 0.37 0.04 0.08]); % behavioral coverage
                    co_ax{10} = axes('units','normalized','Position',[0.35+dir_adj 0.84 0.04 0.08]); %text for n sppikes ,,,
                    %co_ax{11} = axes('units','normalized','Position',[0.12+dir_adj 0.37 0.02 0.08]); % number of co
                    %co_ax{12} = axes('units','normalized','Position',[0.16+dir_adj 0.37 0.02 0.08]); % number of spikes
                    
                    %2D representation:
                    co_ax{11} = axes('units','normalized','Position',[0.05+dir_adj 0.30 0.12 0.17]); % 2D figure allo vs dis
                    co_ax{12} = axes('units','normalized','Position',[0.05+dir_adj 0.47 0.12 0.06]); % ego fr
                    co_ax{17} = axes('units','normalized','Position',[0.17+dir_adj 0.30 0.03 0.17]); % allo fr
                    co_ax{18} = axes('units','normalized','Position',[0.22+dir_adj 0.30 0.03 0.17]); % solo fr
                    co_ax{19} = axes('units','normalized','Position',[0.27+dir_adj 0.30 0.03 0.17]); % solo fr regular bins
                    
                    %                 % stability in time
                    %                 co_ax{13} = axes('units','normalized','Position',[0.05+dir_adj 0.05 0.1 0.12]); % main fig
                    %                 co_ax{14} = axes('units','normalized','Position',[0.05+dir_adj 0.17 0.1 0.04]); % firing rates
                    %                 % information
                    %                 co_ax{15} = axes('units','normalized','Position',[0.2+dir_adj 0.05 0.09 0.12]); % ego information
                    %                 co_ax{16} = axes('units','normalized','Position',[0.32+dir_adj 0.05 0.09 0.12]); % allo information
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% check if dir is signif
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
                        signif=1;
                    end
                    
                    %% info
                    axes(co_ax{10})
                    title(sprintf('# spikes: %d\n# CO: %d',behavior_struct.co(ii_dir).info.n_spikes,behavior_struct.co(ii_dir).info.n_co))
                    box off
                    axis off
                    
                    % main figure
                    axes(co_ax{1})
                    hold on
                    line([0 0],tunnel_limits,'Color',[.7 .7 .7],'LineStyle','--')
                    plot(behavior_struct.co(ii_dir).bsp.dis_m(:),behavior_struct.co(ii_dir).bsp.x_pos(:),'.','color',spike_colors{3})
                    plot(behavior_struct.co(ii_dir).spikes.dis_m(:),behavior_struct.co(ii_dir).spikes.x_pos(:),'.','color',spike_colors{ii_dir},'markersize',8)
                    set(gca,'xlim',[-dis_before_after_co dis_before_after_co],...
                        'ylim',tunnel_limits,...
                        'ytick',tunnel_limits)
                    ylabel('X pos. (m)','fontsize',fsize)
                    xlabel('Dis. between bats (m)','fontsize',fsize)
                    box off
                    
                    % egocentric firing rate
                    axes(co_ax{2})
                    hold on
                    plot(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers,co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate(2:end,:),'color',spike_colors{3},'LineWidth',lwidth);
                    plot(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers,co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate(1,:),'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    max_y = round(max(max(co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate))*1.1*10) / 10;
                    if max_y < 1
                        max_y = 1;
                    end
                    bins_diff = mean(diff(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers));
                    x_limits = [co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(1)-bins_diff/2 co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(end)+bins_diff/2];
                    set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                    str = sprintf('Hz');
                    ylabel(str,'fontsize',fsize)
                    box off
                    
                    % allocentric firing rate
                    axes(co_ax{3})
                    hold on
                    plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                    plot(behavior_struct.co(ii_dir).firing_rate.solo_x_pos(1,:),behavior_struct.co(ii_dir).firing_rate.solo_x_pos(2,:),'color',solo_colors{ii_dir},'LineWidth',2);
                    plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(1,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    max_x = round(max([max(behavior_struct.co(ii_dir).firing_rate.solo_x_pos(1,:)),max(max(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate))])*1.1*10) / 10;
                    if max_x~=0
                        set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                    end
                    str = sprintf('Hz');
                    xlabel(str,'fontsize',fsize)
                    box off
                    
                    % allocentric comparison
                    %                 axes(co_ax{4})
                    %                 hold on
                    %                 bins_centers = behavior_struct.co(ii_dir).solo_co_comparison.bin_centers;
                    %                 n_bins = length(bins_centers);
                    %                 invhibited_ind =behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate  >0;
                    %                 invhibited_bins = zeros(1,n_bins);
                    %                 invhibited_bins(invhibited_ind) =  -behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate(invhibited_ind);
                    %                 excited_bins = zeros(1,n_bins);
                    %                 excited_bins(~invhibited_ind) =  -behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate(~invhibited_ind);
                    %                 barh(bins_centers,invhibited_bins,'BarWidth',1,'facecolor',[0,.8,0])
                    %                 barh(bins_centers,excited_bins,'BarWidth',1,'facecolor',[1,0,0])
                    %                 x_min_max = [min(-behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate) max(-behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate)];
                    %                 if x_min_max(2) == 0
                    %                     x_min_max(2) = 1;
                    %                 end
                    %                 x_limits = [x_min_max(1) - diff(x_min_max)*0.1 , x_min_max(2) + diff(x_min_max)*.3];
                    %                 set(gca,'ycolor',[1 1 1],'ylim',tunnel_limits,'ytick',[],'xlim',x_limits,'xtick',x_min_max,'xticklabel',round(x_min_max*10)/10);
                    %                 adjusted_alpha = alpha / n_bins;
                    %                 sig_p = find(behavior_struct.co(ii_dir).solo_co_comparison.p_diff_firing_rate <=adjusted_alpha);
                    %                 scatter(ones(1,length(sig_p))*x_limits(2),bins_centers(sig_p),[],'k','*')
                    %                 box off
                    %                 str = sprintf('CO-Solo\nHz');
                    %                 xlabel(str,'fontsize',fsize*.7)
                    %                 str = sprintf('n Bins = %d\nAlpha = .05\nadj. Alpha = %.4f',n_bins,adjusted_alpha);
                    %                 [x,y]=ds2nfu([x_limits(2) x_limits(2)],[tunnel_limits(2) tunnel_limits(2)]*1.2);
                    %                 annotation('textarrow',x,y,'string',str,'fontsize',10,'Color',[0,0,0],'HeadStyle','none','LineStyle', 'none', 'TextRotation',0)
                    %
                    % raster of solo spikes
                    axes(co_ax{5})
                    hold on
                    plot(behavior_struct.solo(ii_dir).spikes.ts_usec(:),behavior_struct.solo(ii_dir).spikes.x_pos(:),'.','color','k','markersize',8)
                    time_limits = [min(behavior_struct.solo(ii_dir).bsp.ts_usec(:)) max(behavior_struct.solo(ii_dir).bsp.ts_usec(:))];
                    time_labels=round(round(time_limits*1e-6*(1/60)*10-time_limits(1)*1e-6*(1/60)*10)/10);
                    set(gca,'xlim',time_limits,'xtick',time_limits,'xticklabel',time_labels,'ylim',tunnel_limits,'ytick',[])
                    xlabel ({'Time';'(min)'})
                    title({'Solo';'raster'},'fontsize',10)
                    
                    % solo firing rate - regular bins
                    axes(co_ax{20})
                    hold on
                    %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                    plot(behavior_struct.solo(ii_dir).PSTH_for_field_detection,behavior_struct.solo(1).x_pos_firing_rate{1, 2}   ,'color',solo_colors{ii_dir},'LineWidth',2);
                    plot(behavior_struct.solo(ii_dir).field_height,behavior_struct.solo(ii_dir).field_center,'*r')
                    max_x = 1.1*ceil(max(behavior_struct.solo(ii_dir).PSTH_for_field_detection));
                    if max_x~=0 & ~isnan(max_x)
                        set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                    end
                    str = sprintf('Hz');
                    xlabel(str,'fontsize',fsize)
                    box off
                    
                    
                    
                    % main figure
                    axes(co_ax{11})
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
                        'ylim',ylimits,...
                        'ytick',tunnel_limits)
                    ylabel('X pos. (m)','fontsize',fsize)
                    xlabel('Dis. between bats (m)','fontsize',fsize)
                    max_fr=max(field_density_mat_X_Y(:));
                    text(xlimits(1)+0.4*xlimits(1),0.95*ylimits(2),sprintf('%.1f Hz',max_fr))
                    box off
                    
                    % egocentric firing rate
                    axes(co_ax{12})
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
                    str = sprintf('Hz');
                    ylabel(str,'fontsize',fsize)
                    box off
                    
                    % allocentric firing rate
                    axes(co_ax{17})
                    hold on
                    %  plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                    % plot(behavior_struct.co(ii_dir).firing_rate.solo_x_pos(1,:),behavior_struct.co(ii_dir).firing_rate.solo_x_pos(2,:),'color',solo_colors{ii_dir},'LineWidth',2);
                    plot(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D{bin_dis_i,bin_allo_i},behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers{bin_dis_i,bin_allo_i} ,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    max_x = ceil(max(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D{bin_dis_i,bin_allo_i}) );
                    if max_x~=0
                        set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                    end
                    str = sprintf('Hz');
                    xlabel(str,'fontsize',fsize)
                    box off
                    
                    
                    % solo firing rate -  bins as 2D
                    axes(co_ax{18})
                    hold on
                    %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                    plot(behavior_struct.co(ii_dir).firing_rate.solo_x_pos_same_bin_as_2D{bin_dis_i,bin_allo_i},behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers{bin_dis_i,bin_allo_i},'color',solo_colors{ii_dir},'LineWidth',2);
                    %plot(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D{bin_dis_i,bin_allo_i},behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers{bin_dis_i,bin_allo_i} ,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    max_x = ceil(max(behavior_struct.co(ii_dir).firing_rate.solo_x_pos_same_bin_as_2D{bin_dis_i,bin_allo_i}) );
                    if max_x~=0
                        set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                    end
                    str = sprintf('Hz');
                    xlabel(str,'fontsize',fsize)
                    box off
                    
                    % solo firing rate - regular bins
                    axes(co_ax{19})
                    hold on
                    %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                    plot(behavior_struct.solo(ii_dir).PSTH_for_field_detection,behavior_struct.solo(1).x_pos_firing_rate{1, 2}   ,'color',solo_colors{ii_dir},'LineWidth',2);
                    %sort by field hight-
                    [~,sorted_hights_ind]=sort(behavior_struct.solo(ii_dir).field_height,'descend');
                    fields_to_plot=sorted_hights_ind;
                    num_fileds=4;
                    %take only fields that has enough coverage in the
                    %tuning curve
                    to_skip=[];
                    for fr_i=sorted_hights_ind
                        if sum(isnan(behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1, fr_i}))>=0.5*length(behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1, fr_i})
                            to_skip=[to_skip, fr_i];
                        end
                        
                    end
                    
                    fields_to_plot=setdiff(fields_to_plot,to_skip);
                    if length(fields_to_plot)>num_fileds
                        [~,sorted_hights_ind_ind]=sort(behavior_struct.solo(ii_dir).field_height(fields_to_plot),'descend');
                        fields_to_plot=fields_to_plot(sorted_hights_ind_ind);
                        fields_to_plot=fields_to_plot(1:num_fileds);
                    end
                    if length(behavior_struct.solo(ii_dir).field_height)<=num_fileds
                        num_fileds=length(behavior_struct.solo(ii_dir).field_height);
                        
                    end
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
                    str = sprintf('Hz');
                    xlabel(str,'fontsize',fsize)
                    
                    axis off
                    
                    %% plot per field
                    y_pos_init=[0.05,0.3,0.6];
                    
                    % per field firing rate:
                    if num_fileds>0
                        for smooth_i=1:3
                            sorted_hights_ind_sorted=sort(sorted_hights_ind);
                            dis_y_pos=0.05;
                            
                            count=0;
                            
                            for fr_i=sort(fields_to_plot)
                                count=count+1;
                                r=behavior_struct.co(ii_dir).firing_rate.time_fr_per_field_new_smooth_smooth_window{fr_i}(smooth_i,:);
                                y_pos=y_pos_init(smooth_i)+dis_y_pos*(count-1);
                                axes('units','normalized','Position',[0.35+dir_adj y_pos 0.12 0.05]);
                                if sum(~isnan(r))>=0.5*length(r)
                                    bins_center=time_X_bins_vector_of_centers;
                                    time_signif_field_new_smooth=behavior_struct.co(ii_dir).firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i};
                                    shuf_data=time_signif_field_new_smooth.shuffled_data  ;
                                    
                                    plot(bins_center,shuf_data,'color',spike_colors{3},'LineWidth',lwidth); hold on;
                                    plot(bins_center ,r,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                                    x_limits = [min(bins_center), max(bins_center)];
                                    max_y = round(max(max([r;shuf_data]))*1.1*10) / 10 +1;
                                    if max_y < 1
                                        max_y = 1;
                                    end
                                    if time_signif_field_new_smooth.signif_based_on_extreme_bins==1
                                        neg_signif=time_signif_field_new_smooth.neg_signif;
                                        plot(bins_center(neg_signif),r(neg_signif)+1,'r*')
                                        for width_i=1:length(time_signif_field_new_smooth.width_neg_interp)
                                            if ~isnan(time_signif_field_new_smooth.width_neg_interp(width_i))
                                                plot(time_signif_field_new_smooth.width_line_x_neg_interp(width_i,:),time_signif_field_new_smooth.width_line_y_neg_interp(width_i,:),'k')
                                                text(mean(time_signif_field_new_smooth.width_line_x_neg_interp(width_i,:)),mean(time_signif_field_new_smooth.width_line_y_neg_interp(width_i,:))*0.5,sprintf('%.2f\n%.2f', time_signif_field_new_smooth.width_neg_interp(width_i)./1e6,time_signif_field_new_smooth.neg_rise_time_interp(width_i)./1e6),'color','r')
                                            end
                                        end
                                        pos_signif=time_signif_field_new_smooth.pos_signif;
                                        plot(bins_center(pos_signif),r(pos_signif)+1,'g*')
                                        for width_i=1:length(time_signif_field_new_smooth.width_pos_interp)
                                            if ~isnan(time_signif_field_new_smooth.width_pos_interp(width_i))
                                                plot(time_signif_field_new_smooth.width_line_x_pos_interp(width_i,:),time_signif_field_new_smooth.width_line_y_pos_interp(width_i,:),'k')
                                                text(mean(time_signif_field_new_smooth.width_line_x_pos_interp(width_i,:)),mean(time_signif_field_new_smooth.width_line_y_pos_interp(width_i,:))*0.5,sprintf('%.2f\n%.2f', time_signif_field_new_smooth.width_pos_interp(width_i)./1e6,time_signif_field_new_smooth.pos_rise_time_interp(width_i)./1e6),'color','g')
                                            end
                                        end
                                    end
                                    
                                    if count==1
                                        set(gca,'xlim',x_limits,'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                                    else
                                        set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                                        
                                    end
                                    %                             r=behavior_struct.co(ii_dir).firing_rate.time_x_fr_per_field{1, fr_i};
                                    %                             sparsity=(nanmean(r).^2)/nanmean(r.^2);
                                    %                             std_r=nanstd(r);
                                    %                             mean_r=nanmean(r);
                                    %                             cv=std_r/mean_r;
                                    %                             if behavior_struct.co(ii_dir).firing_rate.time_signif_field{1, fr_i}.signif_based_on_extreme_bins==1
                                    %                                 text(x_limits(2),max_y*0.5,sprintf('SI=%.2f sig bins \n #spikes=%d',behavior_struct.co(ii_dir).firing_rate.time_information_per_spike_per_field{fr_i},behavior_struct.co(ii_dir).firing_rate.number_of_spikes_per_field{fr_i}),'color',[1 0 0])
                                    %                             else
                                    %                                 text(x_limits(2),max_y*0.5,sprintf('SI=%.2f \n #spikes=%d',behavior_struct.co(ii_dir).firing_rate.time_information_per_spike_per_field{fr_i},behavior_struct.co(ii_dir).firing_rate.number_of_spikes_per_field{fr_i}))
                                    %
                                    %                             end
                                    %                             if behavior_struct.co(ii_dir).firing_rate.time_signif_field{1, fr_i}.signif_based_on_CV==1
                                    %
                                    %                                 text(x_limits(2),0,sprintf('cv=%.2f',cv),'color',[1 0 0])
                                    %                             else
                                    %                                 text(x_limits(2),0,sprintf('cv=%.2f',cv))
                                    %                             end
                                    %
                                    % str = sprintf('Hz');
                                    %ylabel(str,'fontsize',fsize)
                                    box off
                                end
                            end
                        end
                        
                    end
                end
                
                
                %% save figure
                
                fig = gcf;
                fig.InvertHardcopy = 'off';
                if ~exist(co_fig_folder_name)
                    mkdir(co_fig_folder_name)
                end
                fig_name=fullfile(co_fig_folder_name,['time_per_field_cell_',num2str(cell_num),'_day_',num2str(day),'_bat_',num2str(bat),'.png']);
                saveas(gcf,fig_name)
                
            end
            
            %% for CO signif cells save also in relevant dir:
            if signif==1
                fig_name=fullfile(co_signif_cells_fig_folder_name,['time_per_field_cell_',num2str(cell_num),'_day_',num2str(day),'_bat_',num2str(bat),'.png']);
                saveas(gcf,fig_name)
            end
            clf(gcf)
            
        else
            continue
        end
    end
    
    %close(gcf)
end
end