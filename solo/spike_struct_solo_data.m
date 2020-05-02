function  solo = spike_struct_solo_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_param_file_name,field_param_file_name)

load(solo_param_file_name)
us_factor=1e6;
%% loop for 2 directions

solo_ind = behavioral_modes.solo_ind;

for ii_dir = 1:2
    
    dir_ind = behavioral_modes.directional_ind{ii_dir};
    
    %find bsp parameters
    bsp_during_solo_ts_usec = bsp_proc_data(tag_i).ts(intersect(solo_ind,dir_ind));
    bsp_during_solo_x_pos = bsp_proc_data(tag_i).pos(intersect(solo_ind,dir_ind),1);
    bsp_during_solo_y_pos = bsp_proc_data(tag_i).pos(intersect(solo_ind,dir_ind),2);
    
    xy=[bsp_during_solo_x_pos,bsp_during_solo_y_pos];
    vel_xy=abs((sqrt(nansum(diff([0 0;xy  ]).^2,2))./diff([0;bsp_during_solo_ts_usec])).*us_factor);
    %find spikes parameters
    epochs=intersect(solo_ind,dir_ind);
    spikes_during_solo_ts_usec = find_spikes_in_epochs (epochs,cell_struct);
    spikes_during_solo_x_pos = interp1(bsp_proc_data(tag_i).ts ,bsp_proc_data(tag_i).pos(:,1),spikes_during_solo_ts_usec);
    spikes_during_solo_y_pos = interp1(bsp_proc_data(tag_i).ts ,bsp_proc_data(tag_i).pos(:,2),spikes_during_solo_ts_usec);
    spikes_during_solo_vel= interp1(bsp_during_solo_ts_usec ,vel_xy,spikes_during_solo_ts_usec);
    
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
     
    data=[bsp_x_pos';bsp_ts_usec';spikes_x_pos';spikes_ts_usec'];
        field_names={'pos','ts','spikes_pos','spikes_ts'};
        FE=cell2struct(data,field_names,1);
       %FE.ts=bsp_x_pos_mat;
        %FE.pos=bsp_x_pos_mat;
        %FE.spikes_ts=spikes_ts_usec_mat;
        %FE.spikes_pos=spikes_x_pos_mat;
        prm.fields=load(field_param_file_name);
    FE_PSTH = FE_compute_PSTH(FE,prm);
%     [timespent_binned, ~, ~, solo_x_pos_firing_rate, ~,~] ...
%         = fn_compute_generic_1D_tuning_new_smooth ...
%         (not_nan_bsp_x_pos, not_nan_spikes_x_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
%    
%     [~, ~, ~, vel_firing_rate, ~,~] ...
%         = fn_compute_generic_1D_tuning_new_smooth ...
%         (vel_xy, spikes_during_solo_vel, vel_bins_vector_of_centers, vel_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
%     
    %calculate velocity curve
     
    %% Field detection
    if ~isempty(not_nan_spikes_x_pos)
        FR_map.PSTH=FE_PSTH.PSTH;
        FR_map.bin_centers=solo_X_bins_vector_of_centers;
        FR_map.bin_size=solo_X_bin_size;
%         data=[bsp_x_pos';bsp_ts_usec';spikes_x_pos';spikes_ts_usec'];
%         field_names={'pos','ts','spikes_pos','spikes_ts'};
%         FE=cell2struct(data,field_names,1);
       %FE.ts=bsp_x_pos_mat;
        %FE.pos=bsp_x_pos_mat;
        %FE.spikes_ts=spikes_ts_usec_mat;
        %FE.spikes_pos=spikes_x_pos_mat;
%         prm.fields=load(field_param_file_name);
        fields=detect_field_based_on_tamir(FR_map,FE,prm) ;
        if ~isempty(fields)
        field_center=[fields.loc];
        field_size=[fields.width_href];
        field_height=[fields.peak];
        field_edges=reshape([fields.edges_prc],2,length([fields.edges_prc])/2)';
        
        PSTH=FR_map.PSTH;
        field_height_norm_by_mean=field_height./nanmean(PSTH);
        else
            field_center=[];
        field_size=[];
        field_height=[];
        field_edges=[];
        field_height_norm_by_mean=[];
        PSTH=FR_map.PSTH;
        end
      
        
 
   %[field_center,field_size,field_height,field_edges,PSTH] = find_fields_PSTH_Basic_withShuffle(flights_struct,solo_param_file_name);
   %field_height_norm_by_mean=field_height./nanmean(PSTH);
    else
        fields=[];
        field_center=[];
        field_size=[];
        field_height=[];
        field_edges=[];
        field_height_norm_by_mean=[];
        PSTH=nan*ones(length(solo_X_bins_vector_of_centers),1);
    end
    %% compute SI:
      information_per_spike = FE_PSTH.SI_bits_spike  ;
    %% save into struct
    solo(ii_dir).bsp.ts_usec = bsp_ts_usec_mat;
    solo(ii_dir).bsp.x_pos = bsp_x_pos_mat;
    solo(ii_dir).bsp.y_pos = bsp_y_pos_mat;
    solo(ii_dir).spikes.ts_usec = spikes_ts_usec_mat;
    solo(ii_dir).spikes.x_pos = spikes_x_pos_mat;
    solo(ii_dir).spikes.y_pos = spikes_y_pos_mat;
    solo(ii_dir).x_pos_firing_rate = {FE_PSTH.PSTH,solo_X_bins_vector_of_centers};
    solo(ii_dir).field_center=field_center;
    solo(ii_dir).field_size=field_size;
    solo(ii_dir).field_height=field_height;
    solo(ii_dir).field_height_norm_by_mean=field_height_norm_by_mean;
    solo(ii_dir).field_edges=field_edges';
    solo(ii_dir).PSTH_for_field_detection=FE_PSTH.PSTH;
    solo(ii_dir).SI=information_per_spike;
    solo(ii_dir).fields=fields;
   % solo(ii_dir).vel_firing_rate=vel_firing_rate;
end

end