function  solo = spike_struct_solo_data (bsp_data,cell_struct,behavioral_modes,tag_i,solo_param_file_name)

load(solo_param_file_name)

%% loop for 2 directions

solo_ind = behavioral_modes.solo_ind;

for ii_dir = 1:2
    
    dir_ind = behavioral_modes.directional_ind{ii_dir};
    
    %find bsp parameters
    bsp_during_solo_ts_usec = bsp_data(tag_i).ts(intersect(solo_ind,dir_ind));
    bsp_during_solo_x_pos = bsp_data(tag_i).pos(intersect(solo_ind,dir_ind),1);
    bsp_during_solo_y_pos = bsp_data(tag_i).pos(intersect(solo_ind,dir_ind),2);
    
    %find spikes parameters
    epochs=intersect(solo_ind,dir_ind);
    spikes_during_solo_ts_usec = find_spikes_in_epochs (epochs,cell_struct);
    spikes_during_solo_x_pos = interp1(bsp_data(tag_i).ts ,bsp_data(tag_i).pos(:,1),spikes_during_solo_ts_usec);
    spikes_during_solo_y_pos = interp1(bsp_data(tag_i).ts ,bsp_data(tag_i).pos(:,2),spikes_during_solo_ts_usec);
    
    % create a matrix from the data, such as each row is one flight
    new_flight_dis_criteria = dis_criteria(ii_dir);
    new_flight_dis_ind = abs(diff(bsp_during_solo_x_pos)) > new_flight_dis_criteria;
    new_flight_time_ind = diff(bsp_during_solo_ts_usec) > new_flight_time_criteria;
    new_flight_ind = find(new_flight_dis_ind & new_flight_time_ind);
    new_flight_start_ind = [1;new_flight_ind+1];
    new_flight_end_ind = [new_flight_ind;length(bsp_during_solo_ts_usec)];
    n_flights = length(new_flight_start_ind);
    
    % prepare empty variables
    bsp_ts_usec = cell(n_flights,1);
    bsp_x_pos = cell(n_flights,1);
    bsp_y_pos = cell(n_flights,1);
    spikes_ts_usec = cell(n_flights,1);
    spikes_x_pos = cell(n_flights,1);
    spikes_y_pos = cell(n_flights,1);
    spikes_ts_usec = cell(n_flights,1);
    bsp_ts_usec= cell(n_flights,1);
    % loop for every flight
    for ii_flight = 1:n_flights
        flight_ind = new_flight_start_ind(ii_flight):new_flight_end_ind(ii_flight);
        bsp_ts_usec{ii_flight} = bsp_during_solo_ts_usec(flight_ind)';
        %bsp_ts_usec{ii_flight}= bsp_during_solo_ts_usec(flight_ind)'.*1e6;
        bsp_x_pos{ii_flight} = bsp_during_solo_x_pos(flight_ind)';
        bsp_y_pos{ii_flight}= bsp_during_solo_y_pos(flight_ind)';
        spikes_ind = (bsp_during_solo_ts_usec(flight_ind(1)) < spikes_during_solo_ts_usec & spikes_during_solo_ts_usec < bsp_during_solo_ts_usec(flight_ind(end)));
        spikes_ts_usec{ii_flight} = spikes_during_solo_ts_usec(spikes_ind);
        spikes_x_pos{ii_flight} = spikes_during_solo_x_pos(spikes_ind);
        spikes_y_pos{ii_flight} = spikes_during_solo_y_pos(spikes_ind);
        %spikes_ts_usec{ii_flight}=spikes_during_solo_ts_usec(spikes_ind).*1e3;
    end
    
    % turn bsp cells into matrices padded with NaNs,
    % each row is a separate flight
    maxSize = max(cellfun(@numel,bsp_ts_usec));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    bsp_ts_usec(cellfun(@isempty,bsp_ts_usec)) = {nan};
    rmat = cellfun(fcn,bsp_ts_usec,'UniformOutput',false);
    bsp_ts_usec_mat = vertcat(rmat{:});
    
    bsp_x_pos(cellfun(@isempty,bsp_x_pos)) = {nan};
    rmat = cellfun(fcn,bsp_x_pos,'UniformOutput',false);
    bsp_x_pos_mat = vertcat(rmat{:});
    
    bsp_y_pos(cellfun(@isempty,bsp_y_pos)) = {nan};
    rmat = cellfun(fcn,bsp_y_pos,'UniformOutput',false);
    bsp_y_pos_mat = vertcat(rmat{:});
    
    % turn spike cells into mat padded with NaNs
    maxSize = max(cellfun(@numel,spikes_ts_usec));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    spikes_ts_usec(cellfun(@isempty,spikes_ts_usec)) = {nan};
    rmat = cellfun(fcn,spikes_ts_usec,'UniformOutput',false);
    spikes_ts_usec_mat = vertcat(rmat{:});
    
    spikes_x_pos(cellfun(@isempty,spikes_x_pos)) = {nan};
    rmat = cellfun(fcn,spikes_x_pos,'UniformOutput',false);
    spikes_x_pos_mat = vertcat(rmat{:});
    
    spikes_y_pos(cellfun(@isempty,spikes_y_pos)) = {nan};
    rmat = cellfun(fcn,spikes_y_pos,'UniformOutput',false);
    spikes_y_pos_mat = vertcat(rmat{:});
    
    % calculate tuning curve
    not_nan_bsp_x_pos = bsp_x_pos_mat(isfinite(bsp_x_pos_mat));
    not_nan_spikes_x_pos = spikes_x_pos_mat(isfinite(spikes_x_pos_mat));
    
    [~, ~, ~, solo_x_pos_firing_rate, ~,~] ...
        = fn_compute_generic_1D_tuning_new_smooth ...
        (not_nan_bsp_x_pos, not_nan_spikes_x_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
    
    %% Field detection
    if ~isempty(not_nan_spikes_x_pos)
        data=[bsp_x_pos';bsp_ts_usec';spikes_x_pos';spikes_ts_usec'];
        field_names={'pos','ts_nlg_usec','spike_pos','spike_ts'};
        flights_struct=cell2struct(data,field_names,1);
        
 
   [field_center,field_size,field_height,field_edges,PSTH] = find_fields_PSTH_Basic_withShuffle(flights_struct,solo_param_file_name);
    else
        field_center=[];
        field_size=[];
        field_height=[];
        field_edges=[];
        PSTH=nan*ones(length(solo_X_bins_vector_of_centers),1);
    end
    
    %% save into struct
    solo(ii_dir).bsp.ts_usec = bsp_ts_usec_mat;
    solo(ii_dir).bsp.x_pos = bsp_x_pos_mat;
    solo(ii_dir).bsp.y_pos = bsp_y_pos_mat;
    solo(ii_dir).spikes.ts_usec = spikes_ts_usec_mat;
    solo(ii_dir).spikes.x_pos = spikes_x_pos_mat;
    solo(ii_dir).spikes.y_pos = spikes_y_pos_mat;
    solo(ii_dir).x_pos_firing_rate = {solo_x_pos_firing_rate,solo_X_bins_vector_of_centers};
    solo(ii_dir).field_center=field_center;
    solo(ii_dir).field_size=field_size;
    solo(ii_dir).field_height=field_height;
    solo(ii_dir).field_edges=field_edges;
    solo(ii_dir).PSTH_for_field_detection=PSTH;
end

end