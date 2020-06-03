function remove_co=find_non_valid_co(CO_point,general_behavior_data_file_name,time_before_after_co_for_co_window,dis_before_after_co,min_dist_opposite_dirs_before_after_CO)
load(general_behavior_data_file_name)
%find length of flight before and after co where the bats are both flying
%in opposite dirs:

%1) Both bats are flying:
both_bats_are_flying=intersect(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
both_bats_are_flying_FE=find(ismember(FE_ind,both_bats_are_flying));
%2) bats with different flight directions:
different_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==-1);

%3) intercet prev conditions:
relevant_flights_FE=intersect(both_bats_are_flying_FE,different_flight_direc_ind_FE);



remove_co=[];
for CO_i=1:length(CO_point)
    
    co_time_usec=ts(CO_point(CO_i));
    bsp_relative_times = ts - co_time_usec;
    bsp_time_criteria_ind = ( - time_before_after_co_for_co_window*us_factor < bsp_relative_times & bsp_relative_times < time_before_after_co_for_co_window*us_factor);
    bsp_dis_criteria_ind = ( abs(bsp_dis_m)<dis_before_after_co)';
    bsp_ind = bsp_time_criteria_ind & bsp_dis_criteria_ind;
    
    
    if ~isempty(relevant_event)
        distance_from_co_before=abs(pos_self_x(distance_change_sign(CO_i))-pos_self_x(different_direc_start(relevant_event)));
        distance_from_co_after=abs(pos_self_x(distance_change_sign(CO_i))-pos_self_x(different_direc_end(relevant_event)));
        if (distance_from_co_before<min_dist_opposite_dirs_before_after_CO) | (distance_from_co_after<min_dist_opposite_dirs_before_after_CO)
            remove_co=[remove_co CO_i];
            
        end
    end
end