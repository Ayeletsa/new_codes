function plot_behav_after_correction(general_behavior_data_file_name,behav_modes_plot,bat,day,co_data)

% %1. plot co that were removed:
% ind=co_data.all_removed_co_ind;
% plot(ts(ind),pos_self_x(ind),'rs')
% 
% %2. plot behavior ind

load(general_behavior_data_file_name)
load(ball_pos_name)
load(behav_params_file_name)

us_factor=1e6;
number_of_plots=9;
%relevant times for plot:
FE_ind=bsp_proc_data(tag_i).flight_ind;
ts=bsp_proc_data(tag_i).ts;

total_time=ts(FE_ind(end))-ts(FE_ind(1));
time_bin=total_time/number_of_plots;
time_vec_for_plot=ts(FE_ind(1)):time_bin:ts(FE_ind(end));

%% figure param:
panel_size=[0.92 0.095];
x_position=0.02;
direc_color=[0 0.447 0.741;0.8500    0.3250    0.0980];

vertical_dist=0.1;
y_position(1)=0.85;
y_position(2)=y_position(1)-vertical_dist;
y_position(3)=y_position(2)-vertical_dist;
y_position(4)=y_position(3)-vertical_dist;
y_position(5)=y_position(4)-vertical_dist;
y_position(6)=y_position(5)-vertical_dist;
y_position(7)=y_position(6)-vertical_dist;
y_position(8)=y_position(7)-vertical_dist;
y_position(9)=y_position(8)-vertical_dist;

%% PLOT
figure('units','normalized','outerposition',[0 0 1 1])

% if tag_i==2
color_self=[0 1 0];
color_other=[1 0 0];

% title:
axes('position', [0.4, y_position(1), 0.1,  0.1]);
title(sprintf('bat %d day %d',bat,day), 'FontSize', 11, 'FontWeight', 'Bold');
box off;
axis off;
% else
%     color_self=[1 0 0];
%     color_other=[0 1 0];
% end
% All behavior
%--------------------------------------------------
behav_color=[0 0 1;0 1 0;1 0 0;0.5 0.5 0.5;1 1 0];

for time_i=1:length(time_vec_for_plot)-1
    ax(time_i)=axes('position',[x_position y_position(time_i) panel_size]);
    relevant_behav_ind=find(ts>=time_vec_for_plot(time_i) & ts<time_vec_for_plot(time_i+1));
   % by_pass=behavioral_modes.bypass_ind(find(ts(behavioral_modes.bypass_ind)>=time_vec_for_plot(time_i) & ts(behavioral_modes.bypass_ind)<time_vec_for_plot(time_i+1)));
    plot(ts(relevant_behav_ind),pos_other_x(relevant_behav_ind),'.','color',color_other)
    hold on;
    plot(ts(relevant_behav_ind),pos_self_x(relevant_behav_ind),'.','color',color_self);
    hold on;
    if sum(isnan(pos_self_x(relevant_behav_ind)))==length(pos_self_x(relevant_behav_ind))
        set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))])
    else
        set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))],'ylim' ,[min(pos_self_x(relevant_behav_ind)) max(pos_self_x(relevant_behav_ind))])
    end
    set(gca,'XTick',[min(ts) max(ts)],'XTickLabel',round(([min(ts) max(ts)])/1e3/60))
 
     %% plot line of dist from ball
  dist_ball_1=ball_1_pos(1)+min_dist_opposite_dirs_before_after_CO+dist_from_the_ball;
  dist_ball_2=ball_2_pos(1)-min_dist_opposite_dirs_before_after_CO-dist_from_the_ball;
  plot(ts([relevant_behav_ind(1) relevant_behav_ind(end)]),[dist_ball_1,dist_ball_1],'k','linewidth',0.5)
  plot(ts([relevant_behav_ind(1) relevant_behav_ind(end)]),[dist_ball_2,dist_ball_2],'k','linewidth',0.5)

  %% plot all behaviors  

 for behav_mod_i=1:length(behav_modes_plot)
        for ii=1:length(behav_modes_plot(behav_mod_i).start)
            start_point=[];
            end_point=[];
            if ts(behav_modes_plot(behav_mod_i).start(ii))>=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))<=time_vec_for_plot(time_i+1)
                start_point=ts(behav_modes_plot(behav_mod_i).start(ii));
                end_point=ts(behav_modes_plot(behav_mod_i).end(ii));
            elseif ts(behav_modes_plot(behav_mod_i).start(ii))>=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))>=time_vec_for_plot(time_i+1)
                start_point=ts(behav_modes_plot(behav_mod_i).start(ii));
                end_point=time_vec_for_plot(time_i+1);
            elseif ts(behav_modes_plot(behav_mod_i).start(ii))<=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))<=time_vec_for_plot(time_i+1)
                start_point=time_vec_for_plot(time_i);
                end_point=ts(behav_modes_plot(behav_mod_i).end(ii));
                
                
            end
            if ~isempty(start_point)
                rec_x=[start_point  end_point  end_point start_point];
                rec_y=[min(pos_self_x)  min(pos_self_x) max(pos_self_x) max(pos_self_x)];
                p=patch(rec_x,rec_y,behav_color(behav_mod_i,:),'EdgeColor','none');
                set(p,'FaceAlpha',0.3)
            end
        end
 end
    %% Plot co there were removed:
    ind=co_data.removed_co_ind_auto;
    relevant_ind_for_time_vec= ind(ts(ind)>=time_vec_for_plot(time_i) & ts(ind)<time_vec_for_plot(time_i+1));
    plot(ts(ind),pos_self_x(ind),'rs')
    %% Plot co there were added:
    
    ind=co_data.manually_added_co;
    relevant_ind_for_time_vec= ind(ts(ind)>=time_vec_for_plot(time_i) & ts(ind)<time_vec_for_plot(time_i+1));
    plot(ts(ind),pos_self_x(ind),'r*')
end
