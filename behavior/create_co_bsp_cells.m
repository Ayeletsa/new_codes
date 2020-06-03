function [remove_co, co_start_end,co_ind,co_bsp_data]=create_co_bsp_data_cells(CO_point,general_behavior_data_file_name,time_before_after_co_for_co_window,dis_before_after_co,min_dist_opposite_dirs_before_after_CO)
load(general_behavior_data_file_name)
%find length of flight before and after co where the bats are both flying
%in opposite dirs:

%1. both bats are flying:
%--------------------------
FE_self=bsp_proc_data(tag_i).flight_ind;
FE_other=bsp_proc_data(3-tag_i).flight_ind;
ind_both_bats_flying=intersect(FE_self,FE_other);
%2. find ind that bats are flying in different directions
%---------------------------------------------------------
temp_bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying,1)';
temp_bsp_other_x_pos = bsp_proc_data(3-tag_i).pos(ind_both_bats_flying,1)';
dir_flight=sign(diff(temp_bsp_x_pos));
dir_flight_other=sign(diff(temp_bsp_other_x_pos));
ind_different_direction=find((dir_flight-dir_flight_other)~=0)+1;
ind_both_bats_flying_different_dirs=ind_both_bats_flying(ind_different_direction);
%% get bsp data valid only for co:
% a. time stamps
bsp_ts = bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs);

% b. x position
bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1)';
other_bat_x_pos=bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1);

% relative distance from the other bat
bsp_dis_m = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1)' - bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1)';

% y position
self_y_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,2)';

% y difference between bats
other_y_pos = bsp_proc_data(3 - tag_i).pos(ind_both_bats_flying_different_dirs,2)';
y_diff = self_y_pos - other_y_pos;

% velocity:
bsp_vel_x=bsp_proc_data(tag_i).vel_x(ind_both_bats_flying_different_dirs)';

bsp_vel_xy=bsp_proc_data(tag_i).vel_xy(ind_both_bats_flying_different_dirs)';

%% cross overs parameters

co_x_positions = bsp_proc_data(tag_i).pos(CO_point,1);
co_times_usec = bsp_proc_data(tag_i).ts(CO_point);
co_directions = sign(co_x_positions - bsp_proc_data(tag_i).pos(CO_point-1,1));
direction_ind = {find(co_directions>0),find(co_directions<0)};

%%
remove_co=[];
co_start_end=[];
co_ind=[];
for ii_dir = 1:2
    
    dir_ind = direction_ind{ii_dir};
    n_co_points = length(dir_ind);
    
    % create empty arrays for CO loop
    bsp_ts_usec = cell(n_co_points,1);
    bsp_time_to_co = cell(n_co_points,1);
    bsp_dis_m_at_co = cell(n_co_points,1);
    bsp_x_pos_at_co = cell(n_co_points,1);
    bsp_y_pos_at_co = cell(n_co_points,1);
    bsp_y_diff_at_co = cell(n_co_points,1);
    bsp_vel_x_at_co = cell(n_co_points,1);
    bsp_vel_xy_at_co = cell(n_co_points,1);

 for ii_co = 1:n_co_points
%for CO_i=1:length(CO_point)
    
    co_time_usec=co_times_usec(dir_ind(ii_co));
    bsp_relative_times = bsp_ts - co_time_usec;
    bsp_time_criteria_ind = ( - time_before_after_co_for_co_window*us_factor < bsp_relative_times & bsp_relative_times < time_before_after_co_for_co_window*us_factor);
    bsp_dis_criteria_ind = ( abs(bsp_dis_m)<dis_before_after_co)';
    bsp_ind = bsp_time_criteria_ind & bsp_dis_criteria_ind;
    relevant_ind=find(bsp_ind);
    
    if sum(bsp_ind)>1
        distance_from_co_before=abs(co_x_positions(dir_ind(ii_co))-bsp_x_pos(relevant_ind(1)));
        distance_from_co_after=abs(co_x_positions(dir_ind(ii_co))-bsp_x_pos(relevant_ind(end)));
        if (distance_from_co_before<min_dist_opposite_dirs_before_after_CO) | (distance_from_co_after<min_dist_opposite_dirs_before_after_CO)
            remove_co=[remove_co dir_ind(ii_co)];
        else
            co_start_end=[co_start_end;relevant_ind(1), relevant_ind(end)];
            co_ind=[co_ind; find(bsp_ind)];
            % b. find samples values
            bsp_ts_usec{ii_co} =  bsp_ts(bsp_ind)';
            bsp_time_to_co{ii_co} = bsp_relative_times(bsp_ind)';
            bsp_dis_m_at_co{ii_co} = bsp_dis_m(bsp_ind);
            bsp_x_pos_at_co{ii_co} = bsp_x_pos(bsp_ind);
            bsp_y_pos_at_co{ii_co} =self_y_pos(bsp_ind);
            bsp_y_diff_at_co{ii_co} = y_diff(bsp_ind);
            bsp_vel_x_at_co{ii_co} = bsp_vel_x(bsp_ind);
            bsp_vel_xy_at_co{ii_co} = bsp_vel_xy(bsp_ind);

        end
    end
 end
co_bsp_data(ii_dir).bsp_ts_usec=bsp_ts_usec;
co_bsp_data(ii_dir).bsp_time_to_co=bsp_time_to_co;
co_bsp_data(ii_dir).bsp_dis_m_at_co=bsp_dis_m_at_co;
co_bsp_data(ii_dir).bsp_x_pos_at_co=bsp_x_pos_at_co;
co_bsp_data(ii_dir).bsp_y_pos_at_co=bsp_y_pos_at_co;
co_bsp_data(ii_dir).bsp_y_diff_co=bsp_y_diff_at_co;
co_bsp_data(ii_dir).bsp_vel_x_co=bsp_vel_x_at_co;
co_bsp_data(ii_dir).bsp_vel_xy_co=bsp_vel_xy_at_co;

end
