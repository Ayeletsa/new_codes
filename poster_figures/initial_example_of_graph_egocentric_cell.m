% plot example for explaining the main graph
clear
cells_num_dir=[296 2];

%% structs' parameters

behavior_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
files = dir(behavior_structs_folder);
behavior_struct_names = {files.name};
co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct\';
files = dir(co_shuffle_structs_folder);
co_shuffle_struct_names = {files.name};


%% figure parameters

spike_colors = {[1  0 1],[1  0 1],[.95 .95 .95]};
solo_colors = {[0 0 0],[0 0 0]};
FR_colors = {[1 1 1],[1 1 1]};%
coherence_colors = {[1 .5 0],[.8 0 .6]};
tunnel_limits = [0 130];
sig_bkg_color = [.8 .95 1];
lwidth = 3;
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

fsize = 14;
y_pos=0.4;
x_solo_psth=0.75;
x_solo_raster=0.695;
x_co_raster=0.5;
x_co_2d=0.7;
x_co_by_place_initial=0.1;


%% plot for each cell

for ii_cell = 1:length(cells_num_dir)
    co_tun_ax = axes('units','normalized','Position',[x_co_raster y_pos+0.175 0.18 0.1]);
    co_raste_ax = axes('units','normalized','Position',[x_co_raster y_pos 0.18 0.17]); %
    co_2d_ax = axes('units','normalized','Position',[x_co_raster y_pos-0.175 0.18 0.17]); %
    solo_raster_ax = axes('units','normalized','Position',[x_solo_raster y_pos 0.05 0.17]); % allo fr comparison
    solo_psth_ax = axes('units','normalized','Position',[x_solo_psth y_pos 0.05 0.17]); % allo fr comparison
    co_by_place_initial_raster_ax = axes('units','normalized','Position',[x_co_by_place_initial y_pos 0.17 0.17]); % allo fr comparison
    co_by_place_initial_psth_ax = axes('units','normalized','Position',[x_co_by_place_initial y_pos+0.175 0.17 0.1]); % allo fr comparison

    %% load data
    cell_file_ind=find(~cellfun(@isempty,regexp(behavior_struct_names,sprintf('cell_%d.mat',cells_num_dir(ii_cell,1)))));
    ii_dir=cells_num_dir(ii_cell,2);
    
    struct_name =behavior_struct_names{cell_file_ind};
    file_name = fullfile(behavior_structs_folder,struct_name);
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
        bin_dis_i=1;
        bin_allo_i=1;
        dis_before_after_co = behavior_struct.co(1).params.dis_before_after_co;
        %% title
        axes('units','normalized','Position',[0.1 0.15 0.7 0.7]);
        text(0.31,0.88,sprintf('Cell %d',cell_num),'fontsize',fsize*1.1,'fontweight','bold')
        text(0.85,0.8, 'Solo data','fontsize',fsize,'fontweight','bold')
        text(0.585,0.8,'Cross-over data','fontsize',fsize,'fontweight','bold')
        text(0.01,0.8,'Cross-over data','fontsize',fsize,'fontweight','bold')

        box off
        axis off
        
        %% co by place initial
        vec_trial=1:size(behavior_struct.co(ii_dir).bsp.dis_m,1);
        mat_trail_bsp=repmat(vec_trial',1,size(behavior_struct.co(ii_dir).bsp.dis_m,2));
       mat_trail_spikes=repmat(vec_trial',1,size(behavior_struct.co(ii_dir).spikes.dis_m,2));

        axes(co_by_place_initial_raster_ax)
        hold on;
        plot(behavior_struct.co(ii_dir).bsp.dis_m,mat_trail_bsp,'.','color',spike_colors{3});
        plot(behavior_struct.co(ii_dir).spikes.dis_m,mat_trail_spikes,'.','color',spike_colors{ii_dir},'markersize',8)

        %plot(cell_co_solo_initial_analysis.co(ii_dir).bsp.x_pos,cell_co_solo_initial_analysis.co(ii_dir).bsp.ts_usec,'.','color',spike_colors{3})
        %plot(cell_co_solo_initial_analysis.co(ii_dir).spikes.x_pos,cell_co_solo_initial_analysis.co(ii_dir).spikes.ts_usec,'.','color',spike_colors{ii_dir},'markersize',8)
        %time_limits = [min(behavior_struct.solo(ii_dir).bsp.ts_usec(:)) max(behavior_struct.solo(ii_dir).bsp.ts_usec(:))];
        %time_limits_str=round([time_limits-time_limits(1)]*1e-6*(1/60));
        trails=size(mat_trail_bsp,1);
        set(gca,'ylim',[0 trails],'ytick',[0 trails])
        %set(gca,'xlim',tunnel_limits,'xtick',tunnel_limits
        set(gca,'xlim',[-dis_before_after_co dis_before_after_co])
        ylabel ('Trails','fontsize',fsize)
        text(0,-15, 'Inter-bat distance (m)','fontsize',fsize,'HorizontalAlignment','center')
         
         %% co by place initial
         axes(co_by_place_initial_psth_ax)
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
        text(-dis_before_after_co-10,max_y/2, str,'fontsize',fsize,'Rotation',90,'HorizontalAlignment','center')
        
        box off
%         hold on;
%           plot(behavior_struct.co(ii_dir).firing_rate.allo_X_bins_vector_of_centers{bin_dis_i,bin_allo_i} ,behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D{bin_dis_i,bin_allo_i},'color',spike_colors{ii_dir},'LineWidth',lwidth);
%                 max_y = ceil(max(behavior_struct.co(ii_dir).firing_rate.allo_x_pos_fr_for_2D{bin_dis_i,bin_allo_i}) );
%                 if max_y~=0
%                     set(gca,'xlim',tunnel_limits,'xTick',[],'ylim',[0 max_y],'ytick',max_y,'color',FR_colors{2});
%                 end
%                 str = sprintf('Hz');
%                 ylabel(str,'fontsize',fsize)
%                 box off
%                 
        %% 2D raster CO
        axes(co_raste_ax)
        hold on
        line([0 0],tunnel_limits,'Color',[.7 .7 .7],'LineStyle','--')
        plot(behavior_struct.co(ii_dir).bsp.dis_m(:),behavior_struct.co(ii_dir).bsp.x_pos(:),'.','color',spike_colors{3})
        plot(behavior_struct.co(ii_dir).spikes.dis_m(:),behavior_struct.co(ii_dir).spikes.x_pos(:),'.','color',spike_colors{ii_dir},'markersize',8)
        set(gca,'xlim',[-dis_before_after_co dis_before_after_co],...
            'ylim',tunnel_limits,...
            'ytick',tunnel_limits(2))
        text(-dis_before_after_co-20,tunnel_limits(2)/2, sprintf('Position\nin tunnel (m)'),'fontsize',fsize,'Rotation',90,'HorizontalAlignment','center')
        set(gca,'xtick',[])
        box off
        %% CO firing rate
        axes(co_tun_ax)
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
        text(-dis_before_after_co-10,max_y/2, str,'fontsize',fsize,'Rotation',90,'HorizontalAlignment','center')
        
        box off
        %% co 2d
        axes(co_2d_ax)
        hold on
        bin_dis_i=1;
        bin_allo_i=1;
       
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
        
        text(-dis_before_after_co-20,tunnel_limits(2)/2, sprintf('Position\nin tunnel (m)'),'fontsize',fsize,'Rotation',90,'HorizontalAlignment','center')
        xlabel('Inter-bat distance (m)','fontsize',fsize)
        max_fr=max(field_density_mat_X_Y(:));
        text(xlimits(2)*0.4,0.95*ylimits(2),sprintf('%.1f Hz',max_fr))
        box off
        
        
        
        %% raster of solo spikes
        axes(solo_raster_ax)
        hold on
        plot(behavior_struct.solo(ii_dir).spikes.ts_usec(:),behavior_struct.solo(ii_dir).spikes.x_pos(:),'.','color','k','markersize',8)
        time_limits = [min(behavior_struct.solo(ii_dir).bsp.ts_usec(:)) max(behavior_struct.solo(ii_dir).bsp.ts_usec(:))];
        time_limits_str=round([time_limits-time_limits(1)]*1e-6*(1/60));
        set(gca,'xlim',time_limits,'xtick',time_limits,'xticklabel',time_limits_str,'ylim',tunnel_limits,'ytick',[])
        xlabel ({'Time';'(min)'},'fontsize',fsize)
        
        %% solo firing rate - regular bins
         axes(solo_psth_ax)
        hold on
        %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
        plot(behavior_struct.solo(ii_dir).PSTH_for_field_detection,behavior_struct.solo(1).x_pos_firing_rate{1, 2},'color',solo_colors{ii_dir},'LineWidth',2);
        max_x = 1.2*ceil(max(behavior_struct.solo(ii_dir).PSTH_for_field_detection));
        plot(behavior_struct.solo(ii_dir).field_height+0.1*max_x,behavior_struct.solo(ii_dir).field_center,'*r')
        if max_x~=0 & ~isnan(max_x)
            set(gca,'ylim',tunnel_limits,'YTick',[],'xlim',[0 max_x],'Xtick',max_x,'color',FR_colors{2});
        end
        str = sprintf('Hz');
        xlabel(str,'fontsize',fsize)
        box off
        
        
        
        
        
        
    end
    
    %% save figure
    
    fig = gcf;
    %fig.InvertHardcopy = 'off';
    dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
    if ~exist(dir_name)
        mkdir(dir_name)
    end
    fig_name=fullfile(dir_name,['initial_egocentric_cell_',num2str(cell_num)]);
    print(gcf, fig_name, '-dpdf',  '-painters');
    
    clf
    
end

