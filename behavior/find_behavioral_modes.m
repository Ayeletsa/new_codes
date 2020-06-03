function [behavioral_modes]=find_behavioral_modes(bsp_data,behav_params_file_name,tag_i,bat,day,behave_ts,dir_param_file_name,ball_pos_name)
%load parameters:
load(behav_params_file_name)
load(dir_param_file_name)

us_factor=1e6;
%% define general behavioral variables:
general_behavior_data_file_name=fullfile(behave_day_struct_folder,['general_behave_analysis_bat',num2str(bat),'_day',day,'.mat']);
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


behav=1;
behav_modes_plot(behav).name='Solo';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(solo_ind,length(pos_self_x));

%% manually remove solo ind:
manually_added_file_name=fullfile(behave_day_struct_folder,['manually_removed_solo_bat_',num2str(bat),'_day_',day,'.mat']);

if correct_manually
    try
        load(manually_added_file_name)
    catch
        %plot:
        CO_point=[];
        removed_co_ind_auto=[];
        correction_type='please remove solo if needed';
        behave_plot(general_behavior_data_file_name,bat,day,CO_point,removed_co_ind_auto,behav_modes_plot,correction_type)
        %correct manually:  
        solo_ind_to_remove=correct_solo(behav_modes_plot(behav).start,behav_modes_plot(behav).end,manual_min_dis_from_CO,ts,pos_self_x);
        save(manually_added_file_name,'solo_ind_to_remove')

    end
    solo_ind(ismember(solo_ind,solo_ind_to_remove))=[];
    [~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(solo_ind,length(pos_self_x));

else
    try
        load(manually_added_file_name)
        solo_ind(solo_ind_to_remove)=[];
        [~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(solo_ind,length(pos_self_x));
        
    catch
        print('no solo correction to load')
    end
end

behavioral_modes.solo_ind=solo_ind;
%% create solo bsp structs:
solo_struct_name=fullfile(behave_solo_struct_folder,['solo_bsp_data_bat_',num2str(bat),'_day_',day]);
[solo_bsp_struct]=create_solo_bsp_data_cells(behavioral_modes,general_behavior_data_file_name,dis_criteria,new_flight_time_criteria,solo_struct_name);
behavioral_modes.solo_bsp_struct=solo_bsp_struct;
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


%% Find potential co events:
[CO_point_auto,removed_co_ind_auto]=find_potential_CO_ind(general_behavior_data_file_name,behav_params_file_name);
co_data.CO_point_auto=CO_point_auto;
co_data.removed_co_ind_auto=removed_co_ind_auto;
%% corect manually co:
manually_added_file_name=fullfile(behave_day_struct_folder,['manually_added_co_bat_',num2str(bat),'_day_',day,'.mat']);
manually_added_co=[];
if correct_manually
    try
        load(manually_added_file_name)
    catch
        % Plot all potential co events:
        correction_type='please add co if needed';
        behave_plot(general_behavior_data_file_name,bat,day,CO_point_auto,removed_co_ind_auto,behav_modes_plot,correction_type)
        % correct manually:
        manually_added_co=mannually_correct_co(pos_self_x,ts,distance_change_sign,manual_min_dis_from_CO);
        save(manually_added_file_name,'manually_added_co')
    end

else
    try
        load(manually_added_file_name)
    catch
        disp('no manually added file')
    end
end

%remove from all
if ~isempty(manually_added_co)
   co_data.manually_added_co=manually_added_co;
    all_point_co_ind=sort([manually_added_co,CO_point_auto]);
else
    co_data.manually_added_co=manually_added_co;
    all_point_co_ind=CO_point_auto;
end

%% find CO window and remove non-valid co
co_struct_name=fullfile(behave_co_struct_folder,['co_bsp_data_bat_',num2str(bat),'_day_',day]);
%remove co that didn't pass threshold of distance
[co_bsp_data,co_ind]=create_co_bsp_data_cells(all_point_co_ind,general_behavior_data_file_name,time_before_after_co_for_co_window,dis_before_after_co,min_dist_opposite_dirs_before_after_CO,co_struct_name,behav_params_file_name,bat,day);
%remove_co_auto=find_non_valid_co(potential_CO_point,general_behavior_data_file_name);

co_data.all_valid_co_ind=all_point_co_ind;
co_data.co_ind=co_ind;
co_data.co_bsp_data=co_bsp_data;
%%
behavioral_modes.co_data=co_data;
behavioral_modes.CO_point=all_point_co_ind;
behavioral_modes.CO_ind=co_ind;

behav=4;
behav_modes_plot(behav).name='CO';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(co_ind,length(pos_self_x));


%% plot final results
plot_behav_after_correction(general_behavior_data_file_name,behav_modes_plot,bat,day,co_data)
file_name=fullfile(behave_analysis_fig_dir_out,['behavioral_modes_bat_',num2str(bat),'_day_',day,'.jpg']);
saveas(gcf,file_name)

%% save all

% %% Cross over:
% %CO_param
% %CO_data
% [CO_ind,CO_point,CO_event]=find_CO_ind(general_behavior_data_file_name,behav_params_file_name);
% 
% %remove empty cells:
% CO_event(cellfun('isempty',CO_event))=[];
% 
% behavioral_modes.CO_event=CO_event;
% behavioral_modes.CO_ind=CO_ind;
% behavioral_modes.CO_point=CO_point;
% 
% % %% plot 1: behavior modes over time
% % manually_corrected_file_name=fullfile(behave_day_struct_folder,['manually_added_behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
% % behavioral_modes=plot_behavior_and_correct_modes(behav_modes_plot,behavioral_modes,bsp_proc_data,tag_i,bat,day,behave_analysis_fig_dir_out,behav_params_file_name,dir_param_file_name,general_behavior_data_file_name,behave_day_struct_folder,manually_corrected_file_name)

%% save struct:

file_name=fullfile(behave_day_struct_folder,['behavioral_modes_bat_',num2str(bat),'_day_',day,'.mat']);
save(file_name,'behavioral_modes')


%% Plot 2: behavioral measurments:
plot_basic_behav_analysis(behavioral_modes,general_behavior_data_file_name,bat,day,behav_params_file_name,behave_analysis_fig_dir_out)
close all
%% CORRECT AND:
% add velocity profile of landing and take offs
% xy of landing and take offs and mark co on them
% velocity profile of co (not just average)
% 





