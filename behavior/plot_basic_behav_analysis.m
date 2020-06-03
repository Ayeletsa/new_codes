function plot_basic_behav_analysis(behavioral_modes,general_behavior_data_file_name,bat,day,behav_params_file_name,behave_analysis_fig_dir_out)
load(general_behavior_data_file_name)
load(behav_params_file_name)
% figure param:
panel_size=[0.1 0.085];
figure('units','normalized','outerposition',[0 0 1 1])
%x pos
horizontal_dist=0.15;
x_position(1)=0.05;
x_position(2)=x_position(1)+horizontal_dist;
x_position(3)=x_position(2)+horizontal_dist;
x_position(4)=x_position(3)+horizontal_dist;
x_position(5)=x_position(4)+horizontal_dist;
x_position(6)=x_position(5)+horizontal_dist;

%y pos
vertical_dist=0.15;
y_position(1)=0.84;
y_position(2)=y_position(1)-vertical_dist;
y_position(3)=y_position(2)-vertical_dist;
y_position(4)=y_position(3)-vertical_dist;
y_position(5)=y_position(4)-vertical_dist;
y_position(6)=y_position(5)-vertical_dist;
y_position(7)=y_position(6)-vertical_dist;
y_position(8)=y_position(7)-vertical_dist;
y_position(9)=y_position(8)-vertical_dist;

% space bins:
X_min=3;
X_max=135;
X_bin_size=5;
X_bins_vector=X_min:X_bin_size:X_max;
X_bins_vector_of_centers=X_bins_vector(1:end-1)+X_bin_size/2;
% space bins CO:
X_bin_size_CO=10;
X_bins_vector_CO=X_min:X_bin_size_CO:X_max;
X_bins_vector_of_centers_CO=X_bins_vector_CO(1:end-1)+X_bin_size_CO/2;

%Time bins:

T_min=ts(FE_ind(1));
T_max=ts(FE_ind(end));
total_time=T_max-T_min;
T_bin_size=total_time/6; %
T_bins_vector=T_min:T_bin_size:T_max;
T_bins_vector_of_centers=T_bins_vector(1:end-1)+T_bin_size/2;

%% plot
% title:
axes('position', [0.4, y_position(1)+0.01, 0.1,  0.1]);
title(sprintf('Bat %d day %d - Total recording time %.1f min, Total flight time %.1f min',bat,day,(total_time/us_factor)/60,(length(FE_ind)/frame_per_second)/60), 'FontSize', 11, 'FontWeight', 'Bold');
box off;
axis off;


%1. Bar graph of time spent in different behaviors
%---------------------------------------------------------------------

% Precentage of flight time
ax1=axes('position',[x_position(1) y_position(1) panel_size]);
behavior_ind=[length(behavioral_modes.solo_ind),length(behavioral_modes.tracking_ind),length(behavioral_modes.CO_ind)];
h=100*(behavior_ind./length(FE_ind));
h(end+1)=100-sum(h);
bar(h)
set(gca,'XTick',1:length(behavior_ind)+1,'XTickLabel',{'Solo','Tracking','CO','None'})
ylabel('% of flight time')
title('Time spent in different behaviors')
xtickangle(ax1,45)
ylim([0 max(h)*1.1])
behavioral_modes.time_spent_in_different_behaviors=h;
% time histograms:
axes('position',[x_position(1) y_position(2) panel_size]);
h=hist(ts(behavioral_modes.solo_ind),T_bins_vector_of_centers);
sharp_ind_solo_hist_time=max(h)./mean(h);
%h=h/frame_per_second;
bar(T_bins_vector_of_centers-ts(FE_ind(1)),h)
set(gca,'XTick',[T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]);
set(gca,'XTickLabel',round([T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]/us_factor/60));
title(sprintf('Solo timings SI=%.2f',sharp_ind_solo_hist_time))
xlabel('Time (minute)')
ylabel('Count')
behavioral_modes.solo_hist_time=h;

axes('position',[x_position(1) y_position(3) panel_size]);
h=hist(ts(behavioral_modes.tracking_ind),T_bins_vector_of_centers);
sharp_ind_tracking_hist_time=max(h)./mean(h);
bar(T_bins_vector_of_centers-ts(FE_ind(1)),h)
set(gca,'XTick',[T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]);
set(gca,'XTickLabel',round([T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]/us_factor/60));
title(sprintf('tracking timings SI=%.2f',sharp_ind_tracking_hist_time))
xlabel('Time (minute)')
ylabel('Count')
behavioral_modes.tracking_hist_time=h;

axes('position',[x_position(1) y_position(5) panel_size]);
h=hist(ts(behavioral_modes.CO_point),T_bins_vector_of_centers);
sharp_ind_CO_hist_time=max(h)./mean(h);
%h=h/frame_per_second;
bar(T_bins_vector_of_centers-ts(FE_ind(1)),h)
set(gca,'XTick',[T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]);
set(gca,'XTickLabel',round([T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]/us_factor/60));
title(sprintf('CO timings SI=%.2f',sharp_ind_CO_hist_time))
xlabel('Time (minute)')
ylabel('Count')
behavioral_modes.CO_hist_time=h;

%%
%2. histogrmas of time spent in different behaviors along the tunnle:
%---------------------------------------------------------------------

%a. solo
axes('position',[x_position(2) y_position(2) panel_size]);
h=hist(pos_self_x(behavioral_modes.solo_ind),X_bins_vector_of_centers);
h=h/frame_per_second;
sharp_ind_solo_hist_space=max(h)./mean(h);
bar(X_bins_vector_of_centers,h)
set(gca,'Xlim',[X_min X_max]);
title(sprintf('Solo position SI=%.2f',sharp_ind_solo_hist_space))
xlabel('X position (m)')
ylabel('Time (s)')
behavioral_modes.solo_hist_space=h;

%b. tracking
axes('position',[x_position(2) y_position(3) panel_size]);
h=hist(pos_self_x(behavioral_modes.tracking_ind),X_bins_vector_of_centers);
h=h/frame_per_second;
sharp_ind_tracking_hist_space=max(h)./mean(h);
bar(X_bins_vector_of_centers,h)
set(gca,'Xlim',[X_min X_max]);
title(sprintf('Tracking position SI=%.2f',sharp_ind_tracking_hist_space))
xlabel('X position (m)')
ylabel('Time (s)')
behavioral_modes.tracking_hist_space=h;


%d. CO
axes('position',[x_position(2) y_position(5) panel_size]);
h=hist(pos_self_x(behavioral_modes.CO_point),X_bins_vector_of_centers_CO);
h_no_edges=h(1+bins_to_remove_from_edge_CO_hist:end-bins_to_remove_from_edge_CO_hist);
sharp_ind_CO_hist_space=max(h_no_edges)./mean(h_no_edges);
bar(X_bins_vector_of_centers_CO,h)
set(gca,'Xlim',[X_min X_max]);
title(sprintf('CO position SI=%.2f',sharp_ind_CO_hist_space))
xlabel('X position (m)')
ylabel('Count')
behavioral_modes.CO_hist_space=h;

%3. velocity XY:
%---------------------------------------------------------------------
%a. velocity profile along the tunnle:
axes('position',[x_position(3) y_position(1) panel_size]);
solo_x=[behavioral_modes.solo_bsp_struct(1).bsp_x_pos;behavioral_modes.solo_bsp_struct(2).bsp_x_pos];    
solo_vel_xy=[behavioral_modes.solo_bsp_struct(1).bsp_vel_xy;behavioral_modes.solo_bsp_struct(2).bsp_vel_xy];
solo_x=cell2mat(solo_x');
solo_vel_xy=cell2mat(solo_vel_xy');

plot(solo_x, solo_vel_xy,'.')
set(gca,'Xlim',[X_min X_max]);
set(gca,'Ylim',[2 9]);
title('Velocity along the tunnle')
xlabel('X position (m)')
ylabel('Velocity (m/s)')
%b. velocity triggered by obstacle

%c. bar of mean velocity during: ALL solo tracked tracking
ax2=axes('position',[x_position(3) y_position(2) panel_size]);
%mean_vel_behavior=[nanmean(abs(velocity_self_xy(FE_ind))),nanmean(abs(velocity_self_xy(solo_ind))),nanmean(abs(velocity_self_xy(tracking_other_ind))),nanmean(abs(velocity_self_xy(being_tracked_ind)))];
mean_vel_behavior=[nanmean(abs(velocity_self_xy(FE_ind))),nanmean(abs(velocity_self_xy(behavioral_modes.solo_ind))),nanmean(abs(velocity_self_xy(behavioral_modes.tracking_ind)))];
[sem_vel(1)]=fn_compute_sem(abs(velocity_self_xy(FE_ind)));
[sem_vel(2)]=fn_compute_sem(abs(velocity_self_xy(behavioral_modes.solo_ind)));
[sem_vel(3)]=fn_compute_sem(abs(velocity_self_xy(behavioral_modes.tracking_ind)));
%[sem(4)]=fn_compute_sem(abs(velocity_self_xy(being_tracked_ind)));
bar(1:length(mean_vel_behavior),mean_vel_behavior); hold on;
errorbar(1:length(mean_vel_behavior),mean_vel_behavior,sem_vel,'.')
set(gca,'XTick',1:3,'XTickLabel',{'All flight','Solo','Tracking'})
xtickangle(ax2,45)
title('velocity (m/s) in differnt behaviors')
behavioral_modes.sem_vel=sem_vel;
behavioral_modes.mean_vel_behavior=mean_vel_behavior;

%d. velocity near the ball
axes('position',[x_position(3) y_position(4) panel_size]);
solo_x=[behavioral_modes.solo_bsp_struct(1).bsp_x_pos;behavioral_modes.solo_bsp_struct(2).bsp_x_pos];    
solo_vel_xy=[behavioral_modes.solo_bsp_struct(1).bsp_vel_xy;behavioral_modes.solo_bsp_struct(2).bsp_vel_xy];
solo_x=cell2mat(solo_x');
solo_vel_xy=cell2mat(solo_vel_xy');
co_x=[behavioral_modes.co_data.co_bsp_data(1).bsp_x_pos_at_co  ;behavioral_modes.co_data.co_bsp_data(2).bsp_x_pos_at_co  ];    
co_vel_xy=[behavioral_modes.co_data.co_bsp_data(1).bsp_vel_xy_at_co  ;behavioral_modes.co_data.co_bsp_data(2).bsp_vel_xy_at_co  ];    
co_x=cell2mat(co_x');
co_vel_xy=cell2mat(co_vel_xy');
load(ball_pos_name)
solo_dist_from_ball=min([abs(solo_x-ball_1_pos(1)); abs(solo_x-ball_2_pos(1))]);
co_dist_from_ball=min([abs(co_x-ball_1_pos(1)); abs(co_x-ball_2_pos(1))]);

plot(solo_dist_from_ball,solo_vel_xy,'k.');
hold on;
plot(co_dist_from_ball,co_vel_xy,'m.');
xlim([0 15])
xlabel('Distance from ball (m)')
ylabel('XY Speed (m/s)')
title('velocity near the balls')

%d. velocity triggered by CO
axes('position',[x_position(3) y_position(5) panel_size]);
all_co_vel=[behavioral_modes.co_data.co_bsp_data(1).bsp_vel_xy_at_co; behavioral_modes.co_data.co_bsp_data(2).bsp_vel_xy_at_co];
all_co_dis=[behavioral_modes.co_data.co_bsp_data(1).bsp_dis_m_at_co; behavioral_modes.co_data.co_bsp_data(2).bsp_dis_m_at_co];
maxSize = max(cellfun(@numel,all_co_vel));
fcn = @(x) [x nan(1,maxSize-numel(x))];
% a. vel
rmat = cellfun(fcn,all_co_vel,'UniformOutput',false);
vel_mat = vertcat(rmat{:});
% b. dis
rmat = cellfun(fcn,all_co_dis,'UniformOutput',false);
dis_mat = vertcat(rmat{:});



plot(dis_mat,vel_mat,'.','color',[.5 .5 .5])
hold on;
plot([0 0],[min(vel_mat(:)) max(vel_mat(:))],'r')
ylim([min(vel_mat(:)) max(vel_mat(:))])
xlim([-40 40])
%ylim([6.5 8.5])
%labels=get(gca,'Xtick');
%set(gca,'Xtick',labels,'XtickLabel',labels)
xlabel('Inter-bat distance (m)')
ylabel('XY Speed (m/s)')
title('velocity triggered by CO')


%4. Y position
%---------------------------------------------------------------------
%a. XY time spent
axes('position',[x_position(4) y_position(1) panel_size]);
y_bin=-1.5:0.25:1.5;
h=hist(pos_self_y,y_bin);
h=100*(h./sum(h));
bar(y_bin,h)
xlim([min(y_bin) max(y_bin)])
ylim([0 max(h)])
xlabel('Y position (m)')
ylabel('% of time')
title('Y time spent')

%b. distance in y during tracking
axes('position',[x_position(4) y_position(3) panel_size]);
y_dist_bin=-2:0.25:2;
h=hist(distnace_other_from_self_y(behavioral_modes.tracking_other_ind),y_dist_bin);
h=100*(h./sum(h));
bar(y_dist_bin,h)
xlim([min(y_dist_bin) max(y_dist_bin)])
xlabel('Y distance (m)')
ylabel('% of time')
title('Distance between the bats during tracking')

%c. distance in y during tracking
axes('position',[x_position(4) y_position(4) panel_size]);
y_dist_bin=-2:0.25:2;
h=hist(distnace_other_from_self_y(behavioral_modes.being_tracked_ind),y_dist_bin);
h=100*(h./sum(h));
bar(y_dist_bin,h)
xlim([min(y_dist_bin) max(y_dist_bin)])
xlabel('Y distance (m)')
ylabel('% of time')
title('Distance between the bats during being tracked')

%e. plot y pos triggered by CO
dir_color=[1,0,0;0,1,0];
axes('position',[x_position(4) y_position(5) panel_size]);
for dir_i=1:2
    co_point_dir=intersect(behavioral_modes.CO_point,behavioral_modes.directional_ind{dir_i});
    mat= repmat([-300:300]',1,length(co_point_dir));
    CO_point_idx = repmat(co_point_dir',[601,1]);
    CO_trigger_idx = CO_point_idx+mat;
    if max(CO_trigger_idx(:))>length(velocity_self_xy)
        CO_trigger_idx(:,end)=[];
    end
    Y_trig_CO=pos_self_y(CO_trigger_idx);
    plot(-300:300,Y_trig_CO,'color',dir_color(dir_i,:))
    hold on;
end
plot([0 0],[min(pos_self_y) max(pos_self_y)],'r')
ylim([-1.7 1.7])
xlim([-300 300])
[x,y]=ds2nfu([150 300],[-1.3 -1.3]);
annotation('arrow',x,y,'color',dir_color(1,:))
[x,y]=ds2nfu([300 150],[-1.5 -1.5]);
annotation('arrow',x,y,'color',dir_color(2,:))
labels=get(gca,'Xtick');
set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
xlabel('Time (s)')
ylabel('y position (m)')
title('Y position triggered by CO')

% plot xy of landing and take offs solo and co
solo_x=[behavioral_modes.solo_bsp_struct(1).bsp_x_pos;behavioral_modes.solo_bsp_struct(2).bsp_x_pos];    
solo_y=[behavioral_modes.solo_bsp_struct(1).bsp_y_pos;behavioral_modes.solo_bsp_struct(2).bsp_y_pos];
solo_x=cell2mat(solo_x');
solo_y=cell2mat(solo_y');
co_x=[behavioral_modes.co_data.co_bsp_data(1).bsp_x_pos_at_co  ;behavioral_modes.co_data.co_bsp_data(2).bsp_x_pos_at_co ];    
co_y=[behavioral_modes.co_data.co_bsp_data(1).bsp_y_pos_at_co  ;behavioral_modes.co_data.co_bsp_data(2).bsp_y_pos_at_co  ];
co_x=cell2mat(co_x');
co_y=cell2mat(co_y');

axes('position',[x_position(5) y_position(1) panel_size]);
plot(solo_x,solo_y,'k.'); hold on;
plot(co_x,co_y,'m.'); hold on;

xlim([min(solo_x)-1, min(solo_x)+10])
xlabel('x position (m)')
ylabel('y position (m)')
title('XY pos near ball 1')

axes('position',[x_position(5) y_position(2) panel_size]);
plot(solo_x,solo_y,'k.'); hold on;
plot(co_x,co_y,'m.'); hold on;
xlim([max(solo_x)-10, max(solo_x)+1])
xlabel('x position (m)')
ylabel('y position (m)')
title('XY pos near ball 2')

% save
fig_name=fullfile(behave_analysis_fig_dir_out,['behavioral_analysis_bat_',num2str(bat),'_day_',num2str(day),'.tif']);
saveas(gcf,fig_name)

clf
