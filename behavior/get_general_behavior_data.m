function get_general_behavior_data(bsp_data,tag_i,general_behavior_data_file_name,ball_pos_name,behave_ts,behav_params_file_name)
us_factor=1e6;
% interpolate position of the other bat to have the same ts for both
% bats
ts=bsp_data(tag_i).ts_us_upsampled;
ts_original=bsp_data(3-tag_i).ts_us_upsampled;
pos_other_original=bsp_data(3-tag_i).pos_upsampled;
pos_other_x=interp1(ts_original,pos_other_original(1,:),ts);
pos_other_y=interp1(ts_original,pos_other_original(2,:),ts);

%define relevant bsp_data struct:
bsp_proc_data(3-tag_i).pos=[pos_other_x;pos_other_y]';
bsp_proc_data(3-tag_i).ts=ts';

bsp_proc_data(tag_i).pos=bsp_data(tag_i).pos_upsampled';
bsp_proc_data(tag_i).ts=bsp_data(tag_i).ts_us_upsampled';

bsp_proc_data=find_flight_ind(bsp_proc_data,behave_ts,behav_params_file_name,ball_pos_name); % find when the bat is flying

FE_ind=bsp_proc_data(tag_i).flight_ind;
ts=bsp_proc_data(tag_i).ts;


%flights position and direction self:
%-----------------------------------------------------------------------------
only_this_bat_is_flying_ind_FE=setdiff(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
only_this_bat_is_flying_ind_FE=find(ismember(FE_ind,only_this_bat_is_flying_ind_FE));
%find the bat X position and velocity
pos_self_x=bsp_proc_data(tag_i).pos(:,1);
pos_self_y=bsp_proc_data(tag_i).pos(:,2);
pos_self_x_FE=pos_self_x(FE_ind);
velocity_self=(diff([0; pos_self_x])./diff([0;bsp_proc_data(tag_i).ts])).*us_factor;
velocity_self(1)=NaN;
velocity_self(isnan(pos_self_x))=NaN;
velocity_self_FE=velocity_self(FE_ind);
% find the xy velocity
pos_self_xy=bsp_proc_data(tag_i).pos(:,1:2);
velocity_self_xy=(sqrt(nansum(diff([0 0;pos_self_xy]).^2,2))./diff([0;bsp_proc_data(tag_i).ts])).*us_factor;
velocity_self_xy(1)=NaN;
velocity_self_xy(isnan(pos_self_x))=NaN;


%flights position and direction other:
%-----------------------------------------------------------------------------
pos_other_x=bsp_proc_data(3-tag_i).pos(:,1);
pos_other_y=bsp_proc_data(3-tag_i).pos(:,2);

pos_other_x_FE=pos_other_x(FE_ind);
velocity_other=(diff([0; bsp_proc_data(3-tag_i).pos(:,1)])./diff([0;bsp_proc_data(3-tag_i).ts])).*us_factor;
velocity_other(1)=NaN;
velocity_other_FE=velocity_other(FE_ind);
velocity_other_FE(only_this_bat_is_flying_ind_FE)=0;
velocity_other(isnan(pos_other_x))=NaN;

%Distance between the bats:
%-----------------------------------------------------------------------------
distnace_other_from_self_FE=pos_other_x_FE-pos_self_x_FE;
distnace_other_from_self=pos_other_x-pos_self_x;
distnace_other_from_self_y=pos_other_y-pos_self_y;

% Both bats are flying:
%-----------------------------------------------------------------------
both_bats_are_flying=intersect(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
both_bats_are_flying_FE=find(ismember(FE_ind,both_bats_are_flying));

%Distance change signs:
%-----------------------------------------------------------------------------
distnace_other_from_self_shifted_plusFE=[0; distnace_other_from_self_FE];
distnace_other_from_self_shifted_minusFE=[distnace_other_from_self_FE; 0];

distance_change_sign=FE_ind(find(sign(distnace_other_from_self_shifted_minusFE.*distnace_other_from_self_shifted_plusFE)==-1)-1);

%save_data
save(general_behavior_data_file_name)

