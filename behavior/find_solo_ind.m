function solo_ind=find_solo_ind(general_behavior_data_file_name,behav_params_file_name)
load(general_behavior_data_file_name)
load(behav_params_file_name)
%Solo is defined when the bats is flying alone or when both are flying but they
%are far away from each other:
%a. only this bat is flying:
a=only_this_bat_is_flying_ind_FE;
% b. both bats are flying but distance is lower than thresh:
b=intersect(both_bats_are_flying_FE, find(abs(distnace_other_from_self_FE)>dist_thresh_solo));
potential_solo_ind=FE_ind(unique([a;b]));

[solo_length,solo_start,solo_end]=find_length_of_consecutive_ind(potential_solo_ind,length(pos_self_x));
long_solo_ind=find(solo_length>min_solo_length);
% take only long flights:
solo_ind=[];
for track_i=1:length(long_solo_ind)
    solo_ind=[solo_ind, solo_start(long_solo_ind(track_i)):solo_end(long_solo_ind(track_i))];
end
