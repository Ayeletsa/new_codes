function bsp_data=find_flight_ind(bsp_data,behave_ts,behav_params_file_name,ball_pos_name);
%find flight index based on distance from ball and velocity:
% flight index is defined:
%1. when the bat is not on the ball(X meters from the ball)
%2. velocity is higher than x
%3. velocity is small but for short times (u-turns for example).

load(behav_params_file_name)
us_factor=1e6;
for tag_i=1:length(bsp_data)
    %load data:
    ts=bsp_data(tag_i).ts;
    pos_x=bsp_data(tag_i).pos(:,1);
    pos_y=bsp_data(tag_i).pos(:,2);
    
    
    %1. find fast flights:
    %-------------------------------------------------------------
    %smooth data before speed compuation
    pos_xy_csaps_pp = csaps(1:length(ts), [pos_x pos_y]', csaps_p);
    pos_xy_csaps = fnval(pos_xy_csaps_pp, 1:length(ts) );
    vel_xy=(sqrt(nansum(diff([0 0;pos_xy_csaps']).^2,2))./diff([0;ts])).*us_factor;
    %remove nan:
    pos_xy_csaps(:,isnan(pos_x)) = nan;
    vel_xy(isnan(pos_x)) = nan;
    
    %remove fast/slow flights:
    vel_xy(vel_xy>flight_upper_vel_thres)=nan;
    fast_flights=find(vel_xy>min_velocity_flight);
    
    %2. Find ind where the bat is not on the ball
    %-------------------------------------------------------------
    if ~exist(ball_pos_name)
        if tag_i==1 %do it just once per day
            figure('units','normalized','outerposition',[0 0 1 1])
            landing_x=pos_x(vel<landing_vel);
            landing_y=pos_y(vel<landing_vel);
            
            % plot the ball area to choose the ball position:
            plot(landing_x,landing_y,'.')
            xlim(ball_1_area);
            title('please choose the ball location')
            ball_1_pos=ginput(1);
            xlim(ball_2_area)
            ball_2_pos=ginput(1);
            
            save(ball_pos_name,'ball_1_pos','ball_2_pos')
            clf;
        end
    else
        load(ball_pos_name)
    end
    ball1_ind=(pos_x<(ball_1_pos(1)+dist_from_the_ball)); % index of measurements on the far ball
    ball2_ind=(pos_x>(ball_2_pos(1)-dist_from_the_ball) ); % index of measurements on the close ball
    not_on_ball=intersect(find(~ball1_ind),find(~ball2_ind));
    
    %4. Take flight only during behavior:
    %-------------------------------------------------------------
    relevant_ind=find(ts>=behave_ts(1) & ts<=behave_ts(2));
    
    %5. unite all conditions:
    %-------------------------------------------------------------
    %relevant_vel_flight=[accepted_slow_flights fast_flights'];
    flight_ind=intersect(fast_flights,not_on_ball);
    flight_ind=intersect(flight_ind,relevant_ind);
    
    
    %6. remove data before and after u-turn:
    %-------------------------------------------------------------
    %find x vel smoothed:
    pos_x_csaps_pp = csaps(1:length(ts), pos_x', csaps_p);
    pos_x_csaps = fnval(pos_x_csaps_pp, 1:length(ts) );
    vel_x=(diff([0; pos_x_csaps'])./diff([0;ts])).*us_factor;
    vel_x(isnan(pos_x)) = nan;
    % find ind of change flight dir:
    velocity_self_shifted_plus_FE=[vel_x(flight_ind); 0];
    velocity_self_shifted_minus_FE=[0 ;vel_x(flight_ind)];
    change_direction_self_ind_FE=find(sign(velocity_self_shifted_plus_FE.*velocity_self_shifted_minus_FE)==-1)-1;
    
    bat_pos_during_ut=pos_x(flight_ind(change_direction_self_ind_FE),1);
    % find distance between ball and co
    dist_to_balls=abs([bat_pos_during_ut-ball_1_pos(1),bat_pos_during_ut-ball_2_pos(1)]);
    % find the minimum of the two distances
    min_dist_to_ball=min(dist_to_balls,[],2);
    UT_potential_ind_FE=change_direction_self_ind_FE(find(min_dist_to_ball>distnace_UT_from_ball));
    ts_FE=bsp_data(tag_i).ts(flight_ind);
    pos_FE=bsp_data(tag_i).pos(flight_ind,:);
    % remove data 1 meter before and after ut:
    FE_ind_to_remove=[];
    for ut_i=1:length(UT_potential_ind_FE)
        ut_time_usec=ts_FE(UT_potential_ind_FE(ut_i));
        ut_pos=pos_FE(UT_potential_ind_FE(ut_i));
        ut_pos_to_remove=[ut_pos-dis_before_after_ut ut_pos+dis_before_after_ut];
        bsp_relative_times = ts_FE - ut_time_usec;
        bsp_time_criteria_ind = ( - time_before_ut_window*us_factor < bsp_relative_times & bsp_relative_times < time_before_ut_window*us_factor);
        bsp_dis_criteria_ind = ( pos_FE(:,1)>ut_pos_to_remove(1) & pos_FE(:,1)<ut_pos_to_remove(2));
        % find if data is around ut:
        if min(ts_FE(bsp_time_criteria_ind))<ut_time_usec & max(ts_FE(bsp_time_criteria_ind))>ut_time_usec
            
            FE_ind_to_remove = [FE_ind_to_remove;find(bsp_time_criteria_ind & bsp_dis_criteria_ind)];
        end
        
    end
    % remove ind before and after ut:
    flight_ind(FE_ind_to_remove)=[];
    
    
    
    %6. save flight ind to struct:
    %-------------------------------------------------------------
    bsp_data(tag_i).flight_ind=flight_ind;
    bsp_data(tag_i).pos_FE=bsp_data(tag_i).pos(flight_ind,:);
    bsp_data(tag_i).ts_FE=bsp_data(tag_i).ts(flight_ind);
    bsp_data(tag_i).vel_x=vel_x;
    bsp_data(tag_i).vel_xy=vel_xy;
    bsp_data(tag_i).vel_x_FE=vel_x(flight_ind);
    bsp_data(tag_i).vel_xy_FE=vel_xy(flight_ind);
    
end