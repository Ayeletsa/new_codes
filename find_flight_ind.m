function bsp_data=find_flight_ind(bsp_data,behave_ts,params_file_name,ball_pos_name);
%find flight index based on distance from ball and velocity:
% flight index is defined:
%1. when the bat is not on the ball(X meters from the ball)
%2. velocity is higher than x
%3. velocity is small but for short times (u-turns for example).

load(params_file_name)
us_factor=1e6;
for tag_i=1:length(bsp_data)
   %load data:
    ts=bsp_data(tag_i).ts;
    pos_x=bsp_data(tag_i).pos(:,1);
    pos_y=bsp_data(tag_i).pos(:,2);
    
    
    %1. find fast flights:
    %-------------------------------------------------------------
    vel=abs((diff([0; bsp_data(tag_i).pos(:,1)])./diff([0;ts])).*us_factor);
    vel(1)=0;
    vel(vel>flight_upper_vel_thres)=nan;
    fast_flights=find(vel>min_velocity_flight);
    
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
    
    %3. find short slow flights:
    %-------------------------------------------------------------
    % we do this because in UT the velocity is low but this is still
    % behavior
    slow_fight_ind=setdiff(1:length(vel),fast_flights);
    %find length of consecutive events:
    [ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(slow_fight_ind,length(vel));
    long_slow_flights=find(ind_length<min_time_for_slow_flights);
    
    accepted_slow_flights=[];
    for flight_i=long_slow_flights
        accepted_slow_flights=[accepted_slow_flights start_ind(flight_i):end_ind(flight_i)];
    end
    
    
    %4. Take flight only during behavior:
    %-------------------------------------------------------------
    relevant_ind=find(ts>=behave_ts(1) & ts<=behave_ts(2));
    
    
    %5. unite all conditions:
    %-------------------------------------------------------------
    relevant_vel_flight=[accepted_slow_flights fast_flights'];
    flight_ind=intersect(relevant_vel_flight,not_on_ball);
    flight_ind=intersect(flight_ind,relevant_ind);
    
    %6. save flight ind to struct:
    %-------------------------------------------------------------
    bsp_data(tag_i).flight_ind=flight_ind;
    bsp_data(tag_i).pos_FE=bsp_data(tag_i).pos(flight_ind,:);
    bsp_data(tag_i).ts_FE=bsp_data(tag_i).ts(flight_ind);
    
end