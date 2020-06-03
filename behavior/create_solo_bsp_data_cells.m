function [solo_bsp_data]=create_solo_bsp_data_cells(behavioral_modes,general_behavior_data_file_name,dis_criteria,new_flight_time_criteria,solo_struct_name)
load(general_behavior_data_file_name)

%% loop for 2 directions

solo_ind = behavioral_modes.solo_ind;

for ii_dir = 1:2
    
    dir_ind = behavioral_modes.directional_ind{ii_dir};
    
    %find bsp parameters
    bsp_during_solo_ts_usec = bsp_proc_data(tag_i).ts(intersect(solo_ind,dir_ind));
    bsp_during_solo_x_pos = bsp_proc_data(tag_i).pos(intersect(solo_ind,dir_ind),1);
    bsp_during_solo_y_pos = bsp_proc_data(tag_i).pos(intersect(solo_ind,dir_ind),2);
    bsp_during_solo_vel_x = bsp_proc_data(tag_i).vel_x(intersect(solo_ind,dir_ind));
    bsp_during_solo_vel_xy = bsp_proc_data(tag_i).vel_xy(intersect(solo_ind,dir_ind));

    % create a matrix from the data, such as each row is one flight
    new_flight_dis_criteria = dis_criteria(ii_dir);
    new_flight_dis_ind = abs(diff(bsp_during_solo_x_pos)) > new_flight_dis_criteria;
    new_flight_time_ind = diff(bsp_during_solo_ts_usec) > new_flight_time_criteria;
    new_flight_ind = find(new_flight_dis_ind | new_flight_time_ind);
    new_flight_start_ind = [1;new_flight_ind+1];
    new_flight_end_ind = [new_flight_ind;length(bsp_during_solo_ts_usec)];
    n_flights = length(new_flight_start_ind);
    
    % prepare empty variables
    bsp_ts_usec = cell(n_flights,1);
    bsp_x_pos = cell(n_flights,1);
    bsp_y_pos = cell(n_flights,1);
    bsp_vel_x= cell(n_flights,1);
    bsp_vel_xy= cell(n_flights,1);
    
    % loop for every flight
    for ii_flight = 1:n_flights
        flight_ind = new_flight_start_ind(ii_flight):new_flight_end_ind(ii_flight);
        bsp_ts_usec{ii_flight} = bsp_during_solo_ts_usec(flight_ind)';
        bsp_x_pos{ii_flight} = bsp_during_solo_x_pos(flight_ind)';
        bsp_y_pos{ii_flight}= bsp_during_solo_y_pos(flight_ind)';
        bsp_vel_x{ii_flight}=bsp_during_solo_vel_x(flight_ind)';
        bsp_vel_xy{ii_flight}=bsp_during_solo_vel_xy(flight_ind)';
    end
    
    solo_bsp_data(ii_dir).bsp_ts_usec=bsp_ts_usec;
    solo_bsp_data(ii_dir).bsp_x_pos=bsp_x_pos;
    solo_bsp_data(ii_dir).bsp_y_pos=bsp_y_pos;
    solo_bsp_data(ii_dir).bsp_vel_x=bsp_vel_x;
    solo_bsp_data(ii_dir).bsp_vel_xy=bsp_vel_xy;
     


solo_struct_name_to_save=[solo_struct_name,'_dir_',num2str(ii_dir),'.mat'];
solo_struct_to_save=solo_bsp_data(ii_dir);
save(solo_struct_name_to_save, '-struct', 'solo_struct_to_save')
end