
% plot imporatnat statistics for cross-overs analysis
close all
clear
clc
cells_num_dir=[290 1; 282 1];







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
tunnel_limits = [0 122];
sig_bkg_color = [.8 .95 1];
lwidth = 2;
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

fsize = 10;

horizontal_dis=0.45;
x_main(1)=0.08;
x_main(2)=x_main(1);

x_behavior=x_main+0.15;
x_spikes=x_behavior+0.06;

vert_dis=0.24;
y_main(1)=0.8;
y_main(2)=y_main(1)-vert_dis;

y_behavior=y_main;
y_spikes=y_main;


for ii_cell = 1:length(cells_num_dir)
              
 
    main_ax{ii_cell} = axes('units','normalized','Position',[x_main(ii_cell) y_main(ii_cell) 0.12 0.12]); %
    behavior_ax{ii_cell} = axes('units','normalized','Position',[x_behavior(ii_cell) y_behavior(ii_cell) 0.03 0.12]);
    
    spikes_ax{ii_cell}=axes('units','normalized','Position',[x_spikes(ii_cell) y_spikes(ii_cell) 0.03 0.12]);
end



%% plot for each cell

for ii_cell = 1:length(cells_num_dir)
    bin_dis_i=1;
    bin_allo_i=1;
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
        
        
        dis_before_after_co = behavior_struct.co(1).params.dis_before_after_co;
        
         % main figure
            axes(main_ax{ii_cell})
            hold on
            line([0 0],tunnel_limits,'Color',[.7 .7 .7],'LineStyle','--')
            plot(behavior_struct.co(ii_dir).bsp.dis_m(:),behavior_struct.co(ii_dir).bsp.x_pos(:),'.','color',spike_colors{3})
            plot(behavior_struct.co(ii_dir).spikes.dis_m(:),behavior_struct.co(ii_dir).spikes.x_pos(:),'.','color',spike_colors{ii_dir},'markersize',8)
            set(gca,'xlim',[-dis_before_after_co dis_before_after_co],...
                'ylim',tunnel_limits,...
                'ytick',tunnel_limits)
                  title(sprintf('Cell %d',cell_num),'fontsize',fsize)

            ylabel('X pos. (m)','fontsize',fsize)
            xlabel('Dis. between bats (m)','fontsize',fsize)
            
            box off
           
        
          % 2D representation
                % bsp only, solo +CO
                axes(behavior_ax{ii_cell})
                hold on
                plot(behavior_struct.co(ii_dir).bsp.y_pos,behavior_struct.co(ii_dir).bsp.x_pos,'.','color',[1 0 1])
                plot(behavior_struct.solo(ii_dir).bsp.y_pos,behavior_struct.solo(ii_dir).bsp.x_pos,'.','color',[0 0 0])
                set(gca,'xlim',[-1.3 1.3],'xtick',[-1.3 1.3],'ylim',tunnel_limits,'ytick',[])
                xlabel({'Y pos.';'(m)'},'fontsize',fsize)
                
                % solo+CO
                axes(spikes_ax{ii_cell})
                hold on
                plot(behavior_struct.solo(ii_dir).bsp.y_pos,behavior_struct.solo(ii_dir).bsp.x_pos,'.','color',[.95 .95 .95])
                plot(behavior_struct.co(ii_dir).bsp.y_pos,behavior_struct.co(ii_dir).bsp.x_pos,'.','color',[1 .9 1])
                plot(behavior_struct.solo(ii_dir).spikes.y_pos,behavior_struct.solo(ii_dir).spikes.x_pos,'.','color',[0 0 0],'markersize',8)
                plot(behavior_struct.co(ii_dir).spikes.y_pos,behavior_struct.co(ii_dir).spikes.x_pos,'.','color',[1 0 1],'markersize',8)
                set(gca,'xlim',[-1.3 1.3],'xtick',[-1.3 1.3],'ylim',tunnel_limits,'ytick',[])
                xlabel({'Y pos.';'(m)'},'fontsize',fsize)
                
              
                
            
        
        
        
    end
end
%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\mid_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,'TwoD_coverage');
print(gcf, fig_name, '-dpdf',  '-painters');

clf


