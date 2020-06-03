function [potential_CO_point,removed_co_ind]=find_potential_CO_ind(general_behavior_data_file_name,behav_params_file_name)
load(general_behavior_data_file_name)
load(behav_params_file_name)

%Cross over is defined when both bats are flying in differnt direction for t time before and t time after
%the distance change sign and te bats are in different directions:

%1) Both bats are flying:
both_bats_are_flying=intersect(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
both_bats_are_flying_FE=find(ismember(FE_ind,both_bats_are_flying));
%2) bats with different flight directions:
different_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==-1);

%3) intercet prev conditions:
relevant_flights_FE=intersect(both_bats_are_flying_FE,different_flight_direc_ind_FE);
%relevant_flights_FE=intersect(relevant_flights_FE,short_distance_ind_FE);

bsp_x_pos = bsp_proc_data(tag_i).pos_FE(relevant_flights_FE,1)';

[different_direc_length,different_direc_start,different_direc_end]=find_length_of_consecutive_ind(FE_ind(relevant_flights_FE),length(pos_self_x));

potential_CO_point=[];
removed_co_ind=[];
co_x_positions = bsp_proc_data(tag_i).pos(distance_change_sign,1);

for CO_i=1:length(distance_change_sign)
    
    relevant_event=find(distance_change_sign(CO_i)>different_direc_start & distance_change_sign(CO_i)<different_direc_end);
    if ~isempty(relevant_event)
        distance_from_co_before=abs(co_x_positions(CO_i)-pos_self_x(different_direc_start(relevant_event)));
        distance_from_co_after=abs(co_x_positions(CO_i)-pos_self_x(different_direc_end(relevant_event)));
        if (distance_from_co_before>min_dist_opposite_dirs_before_after_CO) & (distance_from_co_after>min_dist_opposite_dirs_before_after_CO)
        
            potential_CO_point=[potential_CO_point distance_change_sign(CO_i)];
        else
            removed_co_ind=[removed_co_ind, distance_change_sign(CO_i)];
        end
        
    end
end

