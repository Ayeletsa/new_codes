clear

figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [30 21];
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');
fsize=14;
behav_color=[96 96 96;255 51 255]./255;

x_position=[0.1];
y_position=[0.5 0.2];
panel_size=[0.1 0.1];
%behavior analysis fig for sfn poster
example_day_modes='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\behavioral_modes_bat_2389_day_20180806.mat';
example_day_general='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\general_behave_analysis_bat2389_day20180806.mat';
load(example_day_modes)
load(example_day_general)
load('D:\Ayelet\2bat_proj\Analysis\new_code\params\behav_params.mat')
% space bins:
X_min=3;
X_max=130;
X_bin_size=5;
X_bins_vector=X_min:X_bin_size:X_max;
X_bins_vector_of_centers=X_bins_vector(1:end-1)+X_bin_size/2;
% space bins CO:
X_bin_size_CO=10;
X_bins_vector_CO=X_min:X_bin_size_CO:X_max;
X_bins_vector_of_centers_CO=X_bins_vector_CO(1:end-1)+X_bin_size_CO/2;

%%
axes('position',[x_position y_position(1) panel_size]);
h=hist(pos_self_x(behavioral_modes.solo_ind),X_bins_vector_of_centers);
h=100*h/sum(h);
p=bar(X_bins_vector_of_centers,h,'facecolor',behav_color(1,:));
set(p,'FaceAlpha',0.3)
set(gca,'Xlim',[0 X_max],'Xtick',[0 130]);
set(gca,'ylim',[0 max(h)*1.1]);
title('Solo positions')
xlabel('Position in tunnel (m)')
ylabel('% of time')

axes('position',[x_position y_position(2) panel_size]);
h=hist(pos_self_x(behavioral_modes.CO_point),X_bins_vector_of_centers_CO);
bar(X_bins_vector_of_centers_CO,h,'facecolor',behav_color(2,:))
set(gca,'Xlim',[0 130],'Xtick',[0 130]);
set(gca,'ylim',[0 max(h)*1.1]);
title('Cross-over positions')
xlabel('Position in tunnel (m)')
ylabel('Count')
%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,'behav_example_stat');
print(gcf, fig_name, '-dpdf',  '-painters');

clf
