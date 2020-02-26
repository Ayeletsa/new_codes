
% plot imporatnat statistics for cross-overs analysis
close all
clear
clc
cells_num_dir=[296,2];
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
tunnel_limits = [0 130];
sig_bkg_color = [.8 .95 1];
lwidth = 3;
alpha = .05;


%% open figure

figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [24 21]; % ~A4 size
set(gcf,'DefaultAxesFontSize',9);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');

fsize = 14;

horizontal_dis=0.24;
x_rater(1)=0.08;
x_rater(2)=x_rater(1)+horizontal_dis;
x_rater(3)=x_rater(2)+horizontal_dis;
x_rater(4)=x_rater(3)+horizontal_dis;
x_rater(5)=x_rater(1);
x_rater(6)=x_rater(2);
x_rater(7)=x_rater(3);
x_rater(8)=x_rater(4);
% x_rater(9)=x_rater(1);
% x_rater(10)=x_rater(2);
% x_rater(11)=x_rater(3);
% x_rater(12)=x_rater(4);

x_psth=x_rater;

vert_dis=0.25;
y_raster(1:4)=0.8;
y_raster(5:8)=y_raster(1)-vert_dis;
% y_raster(9:12)=y_raster(5)-vert_dis;


y_psth=y_raster+0.06;

for ii_cell = 1:1
    raster_ax{ii_cell} = axes('units','normalized','Position',[x_rater(ii_cell) y_raster(ii_cell) 0.3 0.06]); % allo fr comparison
    psth_ax{ii_cell} = axes('units','normalized','Position',[x_psth(ii_cell) y_psth(ii_cell) 0.3 0.06]); % allo fr comparison
end




%% plot for each cell

for ii_cell = 1:1
    
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
        
        
        %         %% write cell's ID and stability
        %
        %         dis_before_after_co = behavior_struct.co(1).params.dis_before_after_co;
        %
        %         % write cell's ID
        %         axes('units','normalized','Position',[0.1 0.8 0.18 0.17]);
        %         str = sprintf('Cell #%d, Date: %d, Bat %d',cell_num,day,bat) ;
        %         title(str)
        %         box off;
        %         axis off
        
        
        % raster of solo spikes
        axes(raster_ax{ii_cell})
        hold on
        plot(behavior_struct.solo(ii_dir).spikes.x_pos(:),behavior_struct.solo(ii_dir).spikes.ts_usec(:),'.','color','k','markersize',8)
        time_limits = [min(behavior_struct.solo(ii_dir).bsp.ts_usec(:)) max(behavior_struct.solo(ii_dir).bsp.ts_usec(:))];
        time_limits_str=round([time_limits-time_limits(1)]*1e-6*(1/60));
        x_limits=[8 max(behavior_struct.solo(1).x_pos_firing_rate{1, 2})];
        set(gca,'ylim',time_limits,'ytick',time_limits,'yticklabel',time_limits_str,'xlim',x_limits,'xtick',tunnel_limits)
        str= sprintf('Time\n(min)');
        text(x_limits(1)-0.25*mean(x_limits), time_limits(1)*1.03, str,'fontsize',10,'rotation',90,'VerticalAlignment','bottom')
        xlabel('Position in tunnel (m)')
        
        % solo firing rate - regular bins
        axes(psth_ax{ii_cell})
        hold on
        %plot(co_shuffle_struct(ii_dir).shuffled_data.allo_firing_rate(2:end,:),co_shuffle_struct(ii_dir).shuffled_data.allo_bin_centers,'color',spike_colors{3},'LineWidth',lwidth);
        plot(behavior_struct.solo(1).x_pos_firing_rate{1, 2},behavior_struct.solo(ii_dir).PSTH_for_field_detection   ,'color',solo_colors{ii_dir},'LineWidth',2);
        max_y = 1.2*ceil(max(behavior_struct.solo(ii_dir).PSTH_for_field_detection));
        plot(behavior_struct.solo(ii_dir).field_center,behavior_struct.solo(ii_dir).field_height+0.1*max_y,'*r')
        if max_y~=0 & ~isnan(max_y)
            set(gca,'xlim',x_limits,'XTick',[],'ylim',[0 max_y],'ytick',max_y,'color',FR_colors{2});
        end
        str = sprintf('Hz');
        text(tunnel_limits(1)-0.25*mean(tunnel_limits), max_y/3, str,'fontsize',10,'rotation',90,'VerticalAlignment','bottom')
        box off
        if ii_dir==1
            arrow_point(1,:)=[8 mean(tunnel_limits)/2];
        else
            arrow_point(1,:)=[tunnel_limits(2) (tunnel_limits(2)-mean(tunnel_limits)/2)+8];
            
        end
        arrow_point(2,:)=[max_y*1.1 max_y*1.1];
        [Xf,Yf]=ds2nfu(arrow_point(1,:),arrow_point(2,:));
          annotation('arrow',Xf,Yf)

        title(sprintf('Cell %d',cell_num))
        
    end
end
%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,'place_cells_296');
print(gcf, fig_name, '-dpdf',  '-painters');

clf

    
