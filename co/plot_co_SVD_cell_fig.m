function plot_co_SVD_cell_fig(co_param_file_name,dir_param_file_name,population_param_file_name,per_field_param_file_name)
%% TO DO
% change solo to work on fields struct
% change per field to work on per field struct
% change name of files
% do the same changes also on per field time and inter field codes


%%
% plot imporatnat statistics for cross-overs analysis
load(dir_param_file_name)
load(population_param_file_name)
load(per_field_param_file_name)
load(co_param_file_name)
%% structs' parameters

files = dir(cell_co_solo_initial_analysis_struct_folder);
behavior_struct_names = {files.name};
is_dir=[files.isdir];
behavior_struct_names=behavior_struct_names(~is_dir);

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
for ii_cell =1:length(behavior_struct_names)
    ii_cell
    %try
    signif=0;
    %% load data
    
    struct_name =behavior_struct_names{ii_cell};
    file_name = fullfile(cell_co_solo_initial_analysis_struct_folder,struct_name);
    load(file_name);
    behavior_struct=cell_co_solo_initial_analysis;
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    if isnumeric(day)
        day=num2str(day);
    end
    cell_num=behavior_struct.exp_data.cell_num;
    
    shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d', day,'_c',num2str(cell_num),'.mat'];
    file_name = fullfile(co_shuffle_folder_name,shuffle_struct_name);
    if exist(file_name)
        co_shuffle_struct=load(file_name);
        co_shuffle_struct=co_shuffle_struct.shuffling_struct;
        
        %% write cell's ID and stability
        
        % write cell's ID
        str = sprintf('Cell #%d, Date: %s, Bat %d',cell_num,day,bat) ;
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
        if ~sum(behavior_struct.exp_data.mean_fr)==0
            bar(1,behavior_struct.exp_data.mean_fr  ,.2,'facecolor',[1 .5 0])
            set(gca,'xtick',1,'xticklabel','mean FR','ylim',[0 behavior_struct.exp_data.mean_fr  *1.05],'ytick',behavior_struct.exp_data.mean_fr  ,'yticklabel',round(behavior_struct.exp_data.mean_fr  *10)/10)
        end
        box off
        
        %%     plot data for every direction
        
        for ii_dir = 1:2
            %% tempppppp
            if ~isempty(cell_co_solo_initial_analysis.solo(ii_dir).fields)
                switch per_field_to_plot
                    case 1
                        
                        per_field=cell_co_solo_initial_analysis.co(ii_dir).per_field_href;
                        x=[cell_co_solo_initial_analysis.solo(ii_dir).fields.edges_href];
                        field_edges=reshape(x,2,length(x)/2);
                    case 2
                        per_field=cell_co_solo_initial_analysis.co(ii_dir).per_field_prc;
                        x=[cell_co_solo_initial_analysis.solo(ii_dir).fields.edges_prc];
                        field_edges=reshape(x,2,length(x)/2);
                        
                    case 3
                        per_field=cell_co_solo_initial_analysis.co(ii_dir).per_field_all_spk;
                        x=[cell_co_solo_initial_analysis.solo(ii_dir).fields.edges_all_spk];
                        field_edges=reshape(x,2,length(x)/2);
                end
            else
                field_edges=[];
            end
            %%
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
            
            % info
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
            
            % stability in time
            co_ax{13} = axes('units','normalized','Position',[0.05+dir_adj 0.05 0.1 0.12]); % main fig
            co_ax{14} = axes('units','normalized','Position',[0.05+dir_adj 0.17 0.1 0.04]); % firing rates
            % information
            co_ax{15} = axes('units','normalized','Position',[0.2+dir_adj 0.05 0.09 0.12]); % ego information
            %co_ax{16} = axes('units','normalized','Position',[0.32+dir_adj 0.05 0.09 0.12]); % allo information
            
            % SVD
            co_ax{21} = axes('units','normalized','Position',[0.35+dir_adj 0.58 0.12 0.17]); % all*ego map for SVD
            co_ax{22} = axes('units','normalized','Position',[0.35+dir_adj 0.75 0.12 0.06]); % ego fr
            co_ax{23} = axes('units','normalized','Position',[0.47+dir_adj 0.58 0.03 0.17]); % allo fr
           
            co_ax{24} = axes('units','normalized','Position',[0.35+dir_adj 0.3 0.12 0.17]); % separable map
            
            co_ax{25} = axes('units','normalized','Position',[0.35+dir_adj 0.05 0.12 0.17]); % hist

          
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
            if n_spikes==0
                continue
            else
                % main figure
                axes(co_ax{1})
                hold on
                line([0 0],tunnel_limits,'Color',[.7 .7 .7],'LineStyle','--')
                plot(behavior_struct.co(ii_dir).bsp_full.dis_m(:),behavior_struct.co(ii_dir).bsp_full.x_pos(:),'.','color',spike_colors{3})
                plot(behavior_struct.co(ii_dir).spikes_full.dis_m(:),behavior_struct.co(ii_dir).spikes_full.x_pos(:),'.','color',spike_colors{ii_dir},'markersize',8)
                
                plot(behavior_struct.co(ii_dir).bsp.dis_m(:),behavior_struct.co(ii_dir).bsp.x_pos(:),'.','color',[.5 .6 .6])
                plot(behavior_struct.co(ii_dir).spikes.dis_m(:),behavior_struct.co(ii_dir).spikes.x_pos(:),'.','color','b','markersize',8)
                
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
                plot(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers,behavior_struct.co(ii_dir).firing_rate_full.dis_m(1,:),'color',spike_colors{ii_dir},'LineWidth',lwidth);
                plot(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers,co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate(1,:),'color','b','LineWidth',lwidth);
                
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
                %  plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                plot(behavior_struct.solo(ii_dir).PSTH_for_field_detection,behavior_struct.solo(1).x_pos_firing_rate{1, 2},'color',solo_colors{ii_dir},'LineWidth',2);
                plot(cell_co_solo_initial_analysis.co(ii_dir).firing_rate.allo_x_pos(1,:),cell_co_solo_initial_analysis.co(ii_dir).firing_rate.allo_x_pos (2,:),'color',spike_colors{ii_dir},'LineWidth',lwidth);
                max_x = round(max([max(behavior_struct.solo(ii_dir).PSTH_for_field_detection),max(max(cell_co_solo_initial_analysis.co(ii_dir).firing_rate.allo_x_pos(1,:)))])*1.1*10) / 10;
                if max_x~=0 & ~isnan(max_x)
                    set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                end
                str = sprintf('Hz');
                xlabel(str,'fontsize',fsize)
                box off
                
                % allocentric comparison
                axes(co_ax{4})
                hold on
                if behavior_struct.co(ii_dir).solo_co_comparison.n_bins~=0 & length(behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate)~=1
                    bins_centers = behavior_struct.co(ii_dir).solo_co_comparison.bin_centers;
                    n_bins = length(bins_centers);
                    invhibited_ind =behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate  >0;
                    invhibited_bins = zeros(1,n_bins);
                    invhibited_bins(invhibited_ind) =  -behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate(invhibited_ind);
                    excited_bins = zeros(1,n_bins);
                    excited_bins(~invhibited_ind) =  -behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate(~invhibited_ind);
                    barh(bins_centers,invhibited_bins,'BarWidth',1,'facecolor',[1,0,0])
                    barh(bins_centers,excited_bins,'BarWidth',1,'facecolor',[0,0.8,0])
                    x_min_max = [min(-behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate) max(-behavior_struct.co(ii_dir).solo_co_comparison.diff_mean_firing_rate)];
                    if x_min_max(2) == 0
                        x_min_max(2) = 1;
                    end
                    
                    x_limits = [x_min_max(1) - diff(x_min_max)*0.1 , x_min_max(2) + diff(x_min_max)*.3];
                    set(gca,'ycolor',[1 1 1],'ylim',tunnel_limits,'ytick',[],'xlim',x_limits,'xtick',x_min_max,'xticklabel',round(x_min_max*10)/10);
                    adjusted_alpha = alpha / n_bins;
                    sig_p = find(behavior_struct.co(ii_dir).solo_co_comparison.p_diff_firing_rate <=adjusted_alpha);
                    scatter(ones(1,length(sig_p))*x_limits(2),bins_centers(sig_p),[],'k','*')
                    box off
                    str = sprintf('Solo-CO\nHz');
                    xlabel(str,'fontsize',fsize*.7)
                    str = sprintf('n Bins = %d\nAlpha = .05\nadj. Alpha = %.4f',n_bins,adjusted_alpha);
                    [x,y]=ds2nfu([x_limits(2) x_limits(2)],[tunnel_limits(2) tunnel_limits(2)]*1.2);
                    annotation('textarrow',x,y,'string',str,'fontsize',10,'Color',[0,0,0],'HeadStyle','none','LineStyle', 'none', 'TextRotation',0)
                end
                
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
                fields=behavior_struct.solo(ii_dir).fields;
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
                field_density_mat_X_Y=behavior_struct.co(ii_dir).firing_rate_full.field_density_smoothed_XY_with_NaN;
                
                
                %                 allo_bin_size=behavior_struct.co(ii_dir).firing_rate.allo_bin_size;
                %                 dis_bin_size=behavior_struct.co(ii_dir).firing_rate.dis_bin_size;
                
                [xlimits, ylimits] = fn_plot_2D_field (field_density_mat_X_Y, dis_X_bins_vector_2D, dis_X_bins_vector_of_centers_2D,allo_X_bins_vector_2D, allo_X_bins_vector_of_centers_2D,normalize_to_other_max);
                set(gca,'xlim',xlimits,...
                    'ylim',ylimits,...
                    'ytick',tunnel_limits)
                ylabel('X pos. (m)','fontsize',fsize)
                xlabel('Dis. between bats (m)','fontsize',fsize)
                %x_pos_valid_rectangle=behavior_struct.co(ii_dir).firing_rate.x_pos_valid_rectangle ;
                hold on;
                x_pos_valid=cell_co_solo_initial_analysis.co(ii_dir).x_pos_valid_rectangle ;
                plot(-41*ones(size(x_pos_valid)),x_pos_valid,'.k','markersize',10)
                max_fr=max(field_density_mat_X_Y(:));
                text(xlimits(1)+0.4*xlimits(1),0.95*ylimits(2),sprintf('%.1f Hz',max_fr))
                xlim([-41 40])
                box off
                
                % egocentric firing rate
                axes(co_ax{12})
                hold on
                % plot(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers,co_shuffle_struct(ii_dir).shuffled_data.ego_firing_rate(2:end,:),'color',spike_colors{3},'LineWidth',lwidth);
                plot(dis_X_bins_vector_of_centers_2D  ,behavior_struct.co(ii_dir).firing_rate.dis_x_pos_fr_for_2D,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                max_y = round(max(max(behavior_struct.co(ii_dir).firing_rate.dis_x_pos_fr_for_2D))*1.1*10) / 10;
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
                plot(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D,allo_X_bins_vector_of_centers_2D ,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                max_x = ceil(max(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D) );
                if max_x~=0 &~isnan(max_x)
                    set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                end
                str = sprintf('Hz');
                xlabel(str,'fontsize',fsize)
                box off
                
                
                %
                % solo firing rate - regular bins
                axes(co_ax{19})
                hold on
                %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
                plot(behavior_struct.solo(ii_dir).PSTH_for_field_detection,behavior_struct.solo(1).x_pos_firing_rate{1, 2}   ,'color',solo_colors{ii_dir},'LineWidth',2);
                %sort by field hight-
                [~,sorted_hights_ind]=sort(behavior_struct.solo(ii_dir).field_height,'descend');
                fields_to_plot=1:length(behavior_struct.solo(ii_dir).field_height);
                
                plot(behavior_struct.solo(ii_dir).field_height(fields_to_plot),behavior_struct.solo(ii_dir).field_center(fields_to_plot),'*r')
                if ~isempty(fields_to_plot)
                    for fr_i= fields_to_plot
                        plot([-0.5 -0.5],[field_edges(1,fr_i),field_edges(2,fr_i)],'color',[.5 .5 .5],'linewidth',2)
                    end
                end
                %plot(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D,behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers ,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                max_x = 1.1*ceil(max(behavior_struct.solo(ii_dir).PSTH_for_field_detection));
                if max_x~=0 & ~isnan(max_x)
                    set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[-1 max_x],'Xtick',max_x,'color',FR_colors{2});
                end
                plot([0 0], tunnel_limits,'k')
                str = sprintf('Hz');
                xlabel(str,'fontsize',fsize)
                axis off
                
                
                %% coherence over time
                axes(co_ax{13})
                hold on
                line([0 0],tunnel_limits,'Color',[.7 .7 .7],'LineStyle','--')
                plot(co_shuffle_struct(ii_dir).odd_even_coherence.bsp(2).dis_m(:),co_shuffle_struct(ii_dir).odd_even_coherence.bsp(2).x_pos(:),'.','color',[.9 .9 .9])
                plot(co_shuffle_struct(ii_dir).odd_even_coherence.bsp(1).dis_m(:),co_shuffle_struct(ii_dir).odd_even_coherence.bsp(1).x_pos(:),'.','color',spike_colors{3})
                plot(co_shuffle_struct(ii_dir).odd_even_coherence.spikes(1).dis_m(:),co_shuffle_struct(ii_dir).odd_even_coherence.spikes(1).x_pos(:),'.','color',coherence_colors{1},'markersize',8)
                plot(co_shuffle_struct(ii_dir).odd_even_coherence.spikes(2).dis_m(:),co_shuffle_struct(ii_dir).odd_even_coherence.spikes(2).x_pos(:),'.','color',coherence_colors{2},'markersize',8)
                set(gca,'xlim',[-dis_before_after_co dis_before_after_co],...
                    'ylim',tunnel_limits,...
                    'ytick',tunnel_limits)
                ylabel('X pos. (m)','fontsize',fsize)
                xlabel('Dis. between bats (m)','fontsize',fsize)
                box off
                
                % coherence over time firing rates
                axes(co_ax{14})
                hold on
                plot(co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 3},co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 1},'color',coherence_colors{1},'LineWidth',lwidth);
                plot(co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 3},co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 2},'color',coherence_colors{2},'LineWidth',lwidth);
                max_y = round(max([max(co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 1}),max(co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 2})])*1.1*10) / 10;
                if max_y < 1
                    max_y = 1;
                end
                bins_diff = mean(diff(co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 3}));
                x_limits = [co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 3}(1)-bins_diff/2, co_shuffle_struct(ii_dir).odd_even_coherence.fr{1, 3}(end)+bins_diff/2];
                set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                str = sprintf('Corr. = %.3f', co_shuffle_struct(ii_dir).odd_even_coherence.corr  );
                title(str,'fontsize',10)
                box off
                
                % egocentric information
                axes(co_ax{15})
                if signif_by_SI==1
                    param_name = ('information_per_spike_ego');
                else
                    param_name = ('range');
                end
                hold on
                counts = co_shuffle_struct(ii_dir).shuffled_data.params.(param_name).histogram{1, 1};
                if max(counts) > 0
                    centers = co_shuffle_struct(ii_dir).shuffled_data.params.(param_name).histogram{1, 2};
                    param_value = co_shuffle_struct(ii_dir).shuffled_data.params.(param_name).values(1);
                    bar(centers,counts,1,'FaceColor',[.8 .8 .8],'EdgeColor','none')
                    stairs([centers(1)-(centers(2)-centers(1))/2 centers-(centers(2)-centers(1))/2 centers(length(centers))+(centers(2)-centers(1))/2],[0 counts 0],'LineWidth',1,'color','k')
                    line([param_value param_value],[0 max(counts)],'linestyle','--','color',[1 0 1],'linewidth',1.5);
                    max_x = max([centers(end),param_value]);
                    set(gca,'ylim',[0 max(counts)*1.1],'ytick',max(counts),'yticklabel',round(max(counts)/sum(counts)*100)/100,...
                        'xlim',[0,max_x*1.1],'xtick',[0,max_x*1.1],'xticklabel',round([0,max_x*1.1]*10)/10)
                end
                p_value = co_shuffle_struct(ii_dir).shuffled_data.params.(param_name).p_value;
                if p_value < alpha
                    set(gca, 'color',sig_bkg_color)
                end
                if ~isnan(p_value)
                    str = sprintf('%s = %.3f\n(p = %.4f)\n',param_name,param_value,p_value);
                    title(str,'fontsize',10)
                end
                
                %% SVD
                x=cell_co_solo_initial_analysis.co(ii_dir).firing_rate.field_density_smoothed_XY_with_NaN;
                nan_allo=sum(isnan(x),2)==size(x,2);
                x(nan_allo,:)=[];
                new_allo_vec_of_centers=allo_X_bins_vector_of_centers_2D;
                new_allo_vec_of_centers(nan_allo)=[];
                for_vector=~logical([nan_allo;1]);
                ind_1=find(for_vector);
                for_vector(ind_1(end)+1)=1;
                new_allo_vec=allo_X_bins_vector_2D(for_vector);
                new_tunnel_lim=[min(new_allo_vec) max(new_allo_vec)];
                if isempty(x)
                    continue
                else
                    % guess initial filling of nans:
                    nan_ind=isnan(x);
                    x_no_nan=x;
                    %find nan region:
                    STATS = regionprops(nan_ind,'Centroid','PixelIdxList','ConvexHull');
                    for ii_field=1:1:length(STATS)
                        p = polyshape( STATS(ii_field).ConvexHull); %create a polygon from x y coordinates
                        [cx, cy] = centroid(p); %find the centroid
                        factor = 0.5;
                        %m = scale(p, factor, [cx cy]) %scale p by factor 0.5 w.r.t a point [cx cy]. m's area is 25% of p's area.
                        q = polybuffer(p,+1);
                        polygon_for_plot{ii_field}=q;
                        [xx,yy]=meshgrid(1:size(x,2),1:size(x,1));
                        in = inpolygon(xx(:),yy(:),q.Vertices(:,1),q.Vertices(:,2));
                        ind=sub2ind(size(x),yy(in),xx(in));%FIX!!!!!!!!
                        x_no_nan(STATS(ii_field).PixelIdxList)=nanmean(x_no_nan(ind));
                    end
                    
                    % comute SVD
                    
                    [U,S,V] = svd(x_no_nan);
                    s=diag(S);
                    alpha_svd=1-s(1).^2./sum(s.^2);
                    % shuffle
                    ego_tuning=nanmean(x_no_nan);
                    allo_tuning=nanmean(x_no_nan,2);
                    ego_by_allo_map=ego_tuning.*allo_tuning;
                    ego_by_allo_map=ego_by_allo_map./max(ego_by_allo_map(:));
                    mean_fr=mean(x_no_nan(:));
                    ego_by_allo_map=ego_by_allo_map*mean_fr;
                    for shuffle_i=1:1000
                        r=poissrnd(ego_by_allo_map);
                        gaussian_kernel = fspecial('gaussian',hsize,sigma_a);
                        
                        % Smoothing = convolve with gaussian kernel:
                        r_smoothed = imfilter(r,gaussian_kernel);
                        
                        [U,S,V] = svd(r_smoothed);
                        s=diag(S);
                        alpha_svd_shuffle(shuffle_i)=1-s(1).^2./sum(s.^2);
                    end
                    % plot
                    %----------------------------
                    %ego by allo for SVD
                    axes(co_ax{21})
                    hold on
                    normalize_to_other_max=[];       
                    [xlimits, ylimits] = fn_plot_2D_field (x, dis_X_bins_vector_2D, dis_X_bins_vector_of_centers_2D,new_allo_vec, new_allo_vec_of_centers,normalize_to_other_max);
                    set(gca,'xlim',xlimits,...
                        'ylim',ylimits,...
                        'ytick',ylimits)
                    ylabel('X pos. (m)','fontsize',fsize)
                    xlabel('Dis. between bats (m)','fontsize',fsize)
                    hold on;
                    max_fr=max(field_density_mat_X_Y(:));
                    text(xlimits(1)+0.4*xlimits(1),0.95*ylimits(2),sprintf('%.1f Hz',max_fr))
                    hold on;
%                     for p_i=1:length(polygon_for_plot)
%                         p=polygon_for_plot{p_i};
%                     pgon = polyshape(p.Vertices);
%                     plot(p.Vertices)
%                     end
                    box off
                    % ego tuning
                    axes(co_ax{22})
                    hold on
                    plot(dis_X_bins_vector_of_centers_2D  ,ego_tuning,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    max_y = round(max(max(ego_tuning))*1.1*10) / 10;
                    if max_y < 1
                        max_y = 1;
                    end
                    bins_diff = mean(diff(co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers));
                    x_limits = [co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(1)-bins_diff/2 co_shuffle_struct(ii_dir).shuffled_data.ego_bin_centers(end)+bins_diff/2];
                    set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'Ytick',max_y,'color',FR_colors{1});
                    str = sprintf('Hz');
                    ylabel(str,'fontsize',fsize)
                    box off
                    
                    %allo
                    axes(co_ax{23})
                    hold on
                    plot(allo_tuning,new_allo_vec_of_centers ,'color',spike_colors{ii_dir},'LineWidth',lwidth);
                    max_x = ceil(max(allo_tuning) );
                    if max_x~=0 &~isnan(max_x)
                        set(gca,'ylim',new_tunnel_lim,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
                    end
                    str = sprintf('Hz');
                    xlabel(str,'fontsize',fsize)
                    box off
                    % separable map
                    axes(co_ax{24})
                    hold on
                    normalize_to_other_max=[];
                    [xlimits, ylimits] = fn_plot_2D_field (ego_by_allo_map, dis_X_bins_vector_2D, dis_X_bins_vector_of_centers_2D,new_allo_vec, new_allo_vec_of_centers,normalize_to_other_max);
                    set(gca,'xlim',xlimits,...
                        'ylim',ylimits,...
                        'ytick',ylimits)
                    ylabel('X pos. (m)','fontsize',fsize)
                    xlabel('Dis. between bats (m)','fontsize',fsize)
                    hold on;
                    max_fr=max(field_density_mat_X_Y(:));
                    text(xlimits(1)+0.4*xlimits(1),0.95*ylimits(2),sprintf('%.1f Hz',max_fr))
                    box off
           
                    % hist
                    axes(co_ax{25})
                    hist(alpha_svd_shuffle); hold on; plot([alpha_svd alpha_svd],[0 300],'r')
                    alpha_svd_all=[alpha_svd,alpha_svd_shuffle];
                    alpha_svd_all_sort=sort(alpha_svd_all);
                    p=1-find(alpha_svd_all_sort==alpha_svd)./length(alpha_svd_all_sort);
                    title(sprintf('alpha SVD=%.3f , p=%.3f',alpha_svd,p))
                    
                    
                    
                    
                end
            end
        end
        %% save figure
        
        fig = gcf;
        fig.InvertHardcopy = 'off';
        if ~exist(co_fig_folder_name)
            mkdir(co_fig_folder_name)
        end
        fig_name=fullfile(co_fig_folder_name,['cell_',num2str(cell_num),'_day_',day,'_bat_',num2str(bat),'_signif_by_SI_',num2str(signif_by_SI),'.png']);
        saveas(gcf,fig_name)
        
        %% save fig for n co:
        % less than 10:
        %         n_co_thres=10;
        %         if (behavior_struct.co(1).info.n_co<=n_co_thres) | (behavior_struct.co(2).info.n_co<=n_co_thres)
        %             co_fig_folder_name_thresh=[co_fig_folder_name,'\ten_or_less_co\'];
        %
        %             fig_name=fullfile(co_fig_folder_name_thresh,['cell_',num2str(cell_num),'_day_',day,'_bat_',num2str(bat),'_field_width_type_',num2str(per_field_to_plot),'.png']);
        %             saveas(gcf,fig_name)
        %         end
        %         min_max_thresh=[10 15];
        %         if ((behavior_struct.co(1).info.n_co>min_max_thresh(1))  &  (behavior_struct.co(1).info.n_co<=min_max_thresh(2)) ) | ((behavior_struct.co(2).info.n_co>min_max_thresh(1))  &  (behavior_struct.co(2).info.n_co<=min_max_thresh(2)) )
        %             co_fig_folder_name_thresh=[co_fig_folder_name,'\ten_to_15\'];
        %
        %             fig_name=fullfile(co_fig_folder_name_thresh,['cell_',num2str(cell_num),'_day_',day,'_bat_',num2str(bat),'_field_width_type_',num2str(per_field_to_plot),'.png']);
        %             saveas(gcf,fig_name)
        %         end
        
        %         %% for CO signif cells save also in relevant dir:
        %         if signif==1
        %             fig_name=fullfile(co_signif_cells_fig_folder_name,['cell_',num2str(cell_num),'_day_',num2str(day),'_bat_',num2str(bat),'.png']);
        %             saveas(gcf,fig_name)
        %         end
        clf(gcf)
        
    else
        continue
    end
    %     catch
    %     end
end
close(gcf)

end