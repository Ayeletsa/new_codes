function [UT_potential_ind_FE]=find_UT_ind(general_behavior_data_file_name,behav_params_file_name)


load(general_behavior_data_file_name)
load(behav_params_file_name)


% U-turns are defined bat changed his direction X not near the ball

%1. find ind that bats change directions (vel changed sign):
% %self:
velocity_self_shifted_plus_FE=[velocity_self(FE_ind); 0];
velocity_self_shifted_minus_FE=[0 ;velocity_self(FE_ind)];
change_direction_self_ind_FE=find(sign(velocity_self_shifted_plus_FE.*velocity_self_shifted_minus_FE)==-1)-1;

% UT distance from balls:
load(ball_pos_name)
bat_pos_during_ut=bsp_proc_data(tag_i).pos_FE(change_direction_self_ind_FE,1);
% find distance between ball and co
dist_to_balls=abs([bat_pos_during_ut-ball_1_pos(1),bat_pos_during_ut-ball_2_pos(1)]);
% find the minimum of the two distances
min_dist_to_ball=min(dist_to_balls,[],2);
UT_potential_ind_FE=change_direction_self_ind_FE(find(min_dist_to_ball>distnace_UT_from_ball));

% find events where the diff between potential ut is lower than thresh:
UT_potential_ind_FE(find(diff(UT_potential_ind_FE)<min_diff_UT))=[];


