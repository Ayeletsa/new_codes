function [being_tracked_ind,tracking_other_ind]=find_tracking_ind(general_behavior_data_file_name,behav_params_file_name);
load(general_behavior_data_file_name)
load(behav_params_file_name)


% tracking is defined when the two bats are flying and they are in the same direction (velocity
% with the same sign) and the distance between them is bellow threshold for
% more than t sec

%Conditions:
%1) Both bats are flying:
both_bats_are_flying=intersect(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
both_bats_are_flying_FE=find(ismember(FE_ind,both_bats_are_flying));
%2) short distance:
short_distance_ind_FE=find(abs(distnace_other_from_self_FE)<dist_thresh_tracking);
%3) same flight direction:
same_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==1);

potential_tracking_ind_FE=intersect(same_flight_direc_ind_FE,both_bats_are_flying_FE);
potential_tracking_ind_FE=intersect(potential_tracking_ind_FE,short_distance_ind_FE);

%4) divide to when the bat is tracking of being tracked:
%a. relative distance between the bats:
bat_behind=FE_ind([find(distnace_other_from_self_FE>0 & velocity_self_FE>0);find(distnace_other_from_self_FE<0 & velocity_self_FE<0)]);
bat_ahead=FE_ind([find(distnace_other_from_self_FE<0 & velocity_self_FE>0);find(distnace_other_from_self_FE>0 & velocity_self_FE<0)]);
%b. combine with previous conditions:
being_tracked_ind_also_short=intersect(FE_ind(potential_tracking_ind_FE),bat_ahead);
tracking_other_ind_also_short=intersect(FE_ind(potential_tracking_ind_FE),bat_behind);

%5) minimum time for tracking:
%initialize to find start and end of tracking:
%a. being tracked:
%-------------------------------
[being_tracked_length,being_tracked_start,being_tracked_end]=find_length_of_consecutive_ind(being_tracked_ind_also_short,length(pos_self_x));
long_being_tracked_ind=find(being_tracked_length>min_tracking_length);
being_tracked_ind=[];
for track_i=1:length(long_being_tracked_ind)
    being_tracked_ind=[being_tracked_ind, being_tracked_start(long_being_tracked_ind(track_i)):being_tracked_end(long_being_tracked_ind(track_i))];
end

%b. tracking
%-------------------------------
[tracking_other_length,tracking_other_start,tracking_other_end]=find_length_of_consecutive_ind(tracking_other_ind_also_short,length(pos_self_x));
long_tracking_other_ind=find(tracking_other_length>min_tracking_length);
tracking_other_ind=[];
for track_i=1:length(long_tracking_other_ind)
    tracking_other_ind=[tracking_other_ind, tracking_other_start(long_tracking_other_ind(track_i)):tracking_other_end(long_tracking_other_ind(track_i))];
end
