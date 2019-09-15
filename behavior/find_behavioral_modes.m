function [behavioral_modes]=find_behavioral_modes(bsp_data,behav_params_file_name,tag_i,bat,day,behave_ts,dir_param_file_name,ball_pos_name)
%load parameters:
load(behav_params_file_name)
load(dir_param_file_name)

us_factor=1e6;
%% define general behavioral variables:
general_behavior_data_file_name=fullfile(behave_day_struct_folder,['general_behave_analysis_bat',num2str(bat),'_day',num2str(day),'.mat']);
get_general_behavior_data(bsp_data,tag_i,general_behavior_data_file_name,ball_pos_name,behave_ts,behav_params_file_name)
 
load(general_behavior_data_file_name)

behavioral_modes.pos_other=bsp_data(3-tag_i).pos;
behavioral_modes.pos_self=bsp_data(tag_i).pos;
behavioral_modes.ts=ts;
% find flight directions:
behavioral_modes.directional_ind{1}=FE_ind(find(velocity_self_FE>0));
behavioral_modes.directional_ind{2}=FE_ind(find(velocity_self_FE<0));


%% Solo

solo_ind=find_solo_ind(general_behavior_data_file_name,behav_params_file_name);

behavioral_modes.solo_ind=solo_ind;

behav=1;
behav_modes_plot(behav).name='Solo';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(solo_ind,length(pos_self_x));

%% Tracking
[being_tracked_ind,tracking_other_ind]=find_tracking_ind(general_behavior_data_file_name,behav_params_file_name);

% assign to structs:
behavioral_modes.being_tracked_ind=being_tracked_ind;
behavioral_modes.tracking_other_ind=tracking_other_ind;
%combine:
behavioral_modes.tracking_ind=sort([tracking_other_ind,being_tracked_ind]);
tracking_ind=behavioral_modes.tracking_ind;

% Find ind of bypass-
%-------------------------------------------------------------------------
bypass_ind=find_bypass_ind(distnace_other_from_self,tracking_ind);

behavioral_modes.bypass_ind=bypass_ind;

%find start_end for plots:
behav=2;
behav_modes_plot(behav).name='tracking';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(behavioral_modes.tracking_other_ind,length(pos_self_x));
behav=3;
behav_modes_plot(behav).name='being tracked';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(behavioral_modes.being_tracked_ind,length(pos_self_x));

%% Cross over:
%CO_param
%CO_data
[CO_ind,CO_point,CO_event]=find_CO_ind(general_behavior_data_file_name,behav_params_file_name);

%remove empty cells:
CO_event(cellfun('isempty',CO_event))=[];

behavioral_modes.CO_event=CO_event;
behavioral_modes.CO_ind=CO_ind;
behavioral_modes.CO_point=CO_point;
%% U-turns
[CO_no_UT,UT_ind]=find_UT_ind(general_behavior_data_file_name,behav_params_file_name,CO_ind);

behavioral_modes.CO_no_UT=CO_no_UT;
behavioral_modes.UT_ind=UT_ind;

behav=4;
behav_modes_plot(behav).name='CO';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(CO_no_UT,length(pos_self_x));
behav=5;
behav_modes_plot(behav).name='UT';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(UT_ind,length(pos_self_x));

%% plot 1: behavior modes over time
plot_behavior_and_correct_modes(behav_modes_plot,behavioral_modes,bsp_proc_data,tag_i,bat,day,behave_analysis_fig_dir_out,behav_params_file_name,dir_param_file_name,general_behavior_data_file_name,behave_day_struct_folder)

%% save struct:

file_name=fullfile(behave_day_struct_folder,['behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
save(file_name,'behavioral_modes')


%% Plot 2: behavioral measurments:
plot_basic_behav_analysis(behavioral_modes,general_behavior_data_file_name,bat,day,behav_params_file_name,behave_analysis_fig_dir_out)






