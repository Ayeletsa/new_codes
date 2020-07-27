
function data=PV_get_y_intersect_data_velocity(data,vel_reduction_thresh)
% exclude data where the velocity eas slow lower than threshold
for dir_i=1:2
    ind_intersect=struct();
    other_dir_ind_intersect=struct();
    
    % find reduction in velocity in co in short distance
    vel_prop_to_median=data(dir_i).co.bsp.vel_prop_to_median;
    valid_flights_based_on_vel=vel_prop_to_median>=vel_reduction_thresh;
    %
    
    %keep only valid CO based on velocity
    %bsp:
    data(dir_i).co.bsp.flight_ts=data(dir_i).co.bsp.flight_ts(valid_flights_based_on_vel,:);
    data(dir_i).co.bsp.flight_x_pos=data(dir_i).co.bsp.flight_x_pos(valid_flights_based_on_vel,:);
    data(dir_i).co.bsp.flight_y_pos=data(dir_i).co.bsp.flight_y_pos(valid_flights_based_on_vel,:);
    data(dir_i).co.bsp.flight_dis=data(dir_i).co.bsp.flight_dis(valid_flights_based_on_vel,:);
    
    %spikes:
    data(dir_i).co.spikes.flight_ts=data(dir_i).co.spikes.flight_ts(valid_flights_based_on_vel,:);
    data(dir_i).co.spikes.flight_x_pos=data(dir_i).co.spikes.flight_x_pos(valid_flights_based_on_vel,:);
    data(dir_i).co.spikes.flight_y_pos=data(dir_i).co.spikes.flight_y_pos(valid_flights_based_on_vel,:);
    data(dir_i).co.spikes.flight_dis=data(dir_i).co.spikes.flight_dis(valid_flights_based_on_vel,:);
    
    
    
end




