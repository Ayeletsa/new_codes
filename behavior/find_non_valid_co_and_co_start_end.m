function [remove_co, co_start_end,co_ind]=find_non_valid_co_and_co_start_end(CO_point,general_behavior_data_file_name,time_before_after_co_for_co_window,dis_before_after_co,min_dist_opposite_dirs_before_after_CO)
load(general_behavior_data_file_name)
%find length of flight before and after co where the bats are both flying
%in opposite dirs:

%1) Both bats are flying:
both_bats_are_flying=intersect(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
%2) bats with different flight directions:
different_flight_direc_ind=find(sign(velocity_other.*velocity_self)==-1);

%3) intersect prev conditions:
relevant_flight_ind=intersect(both_bats_are_flying,different_flight_direc_ind);

logical_flight_ind=zeros(1,length(ts));
logical_flight_ind(relevant_flight_ind)=1;
relevant_ts=ts;
relevant_ts(~logical_flight_ind)=nan;
relevant_pos_self_x=pos_self_x;
relevant_pos_self_x(~logical_flight_ind)=nan;
releval_distnace_other_from_self=distnace_other_from_self;
releval_distnace_other_from_self(~logical_flight_ind)=nan;
relevant_pos_self_y=pos_self_y;
relevant_pos_self_y(~logical_flight_ind)=nan;
relevant_pos_self_y=pos_self_y;
relevant_pos_self_y(~logical_flight_ind)=nan;


remove_co=[];
co_start_end=[];
co_ind=[];
for CO_i=1:length(CO_point)
    
    co_time_usec=ts(CO_point(CO_i));
    bsp_relative_times = relevant_ts - co_time_usec;
    bsp_time_criteria_ind = ( - time_before_after_co_for_co_window*us_factor < bsp_relative_times & bsp_relative_times < time_before_after_co_for_co_window*us_factor);
    bsp_dis_criteria_ind = ( abs(releval_distnace_other_from_self)<dis_before_after_co)';
    bsp_ind = bsp_time_criteria_ind & bsp_dis_criteria_ind';
    relevant_ind=find(bsp_ind);
    
    if sum(bsp_ind)>1
        distance_from_co_before=abs(relevant_pos_self_x(CO_point(CO_i))-relevant_pos_self_x(relevant_ind(1)));
        distance_from_co_after=abs(relevant_pos_self_x(CO_point(CO_i))-relevant_pos_self_x(relevant_ind(end)));
        if (distance_from_co_before<min_dist_opposite_dirs_before_after_CO) | (distance_from_co_after<min_dist_opposite_dirs_before_after_CO)
            remove_co=[remove_co CO_i];
        else
            co_start_end=[co_start_end;relevant_ind(1), relevant_ind(end)];
            co_ind=[co_ind; find(bsp_ind)];
            
        end
    end
end