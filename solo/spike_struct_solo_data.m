function  solo = spike_struct_solo_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_param_file_name,field_param_file_name,solo_struct_name)

load(solo_param_file_name)
us_factor=1e6;

%load spikes:
spikes_ts = cell_struct.spikes.spikes_ts_usec;

%find spike ts that are not during bsp nan:
spike_pos=interp1(bsp_proc_data(tag_i).ts,bsp_proc_data(tag_i).pos(:,1),spikes_ts);
spikes_ts(isnan(spike_pos))=[];
%% loop for 2 directions

for ii_dir = 1:2
    %load behavioral data for solo:
    solo_struct_name_to_load=[solo_struct_name,'_dir_',num2str(ii_dir),'.mat'];
    load(solo_struct_name_to_load)
    %% remove flights for unstabel cells:
    [logical_vec_of_active_flight]=correct_behavior_to_work_only_for_active_time_per_cell(cell_struct,bsp_ts_usec);
    ind_of_relevant_flight=find(logical_vec_of_active_flight);
    bsp_ts_usec=bsp_ts_usec(ind_of_relevant_flight);
    bsp_x_pos=bsp_x_pos(ind_of_relevant_flight);
    bsp_y_pos=bsp_y_pos(ind_of_relevant_flight);
    bsp_vel_x=bsp_vel_x(ind_of_relevant_flight);
    bsp_vel_xy=bsp_vel_xy(ind_of_relevant_flight);
    
    solo(ii_dir).bsp.ind_of_relevant_flight = ind_of_relevant_flight;

    %%
    start_ts=mat2cell([cellfun(@(v) v(1),bsp_ts_usec)],ones(length(bsp_ts_usec),1));
    end_ts=mat2cell([cellfun(@(v) v(end),bsp_ts_usec)],ones(length(bsp_ts_usec),1));
%% loop over flights to find spikes per flight:
    n_flights=length(bsp_ts_usec);

    spikes_x_pos = cell(n_flights,1);
    spikes_y_pos = cell(n_flights,1);
    spikes_ts_usec = cell(n_flights,1);
    % loop for every flight
    for ii_flight = 1:n_flights
        spikes_ind = find(spikes_ts> bsp_ts_usec{ii_flight}(1) & spikes_ts<bsp_ts_usec{ii_flight}(end));

        spikes_ts_usec{ii_flight} = spikes_ts(spikes_ind);
        spikes_x_pos{ii_flight} = interp1(bsp_ts_usec{ii_flight},bsp_x_pos{ii_flight},spikes_ts_usec{ii_flight});
        spikes_y_pos{ii_flight} = interp1(bsp_ts_usec{ii_flight},bsp_y_pos{ii_flight},spikes_ts_usec{ii_flight});
       
    end
    %% convert cells to mats:
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
    
    %% calculate tuning curve
    not_nan_bsp_x_pos = bsp_x_pos_mat(isfinite(bsp_x_pos_mat));
    not_nan_spikes_x_pos = spikes_x_pos_mat(isfinite(spikes_x_pos_mat));
    %create data for computing PSTH:     
    data=[bsp_x_pos';bsp_ts_usec';spikes_x_pos';spikes_ts_usec';start_ts';end_ts'];
    field_names={'pos','ts','spikes_pos','spikes_ts','start_ts','end_ts'};
    FE=cell2struct(data,field_names,1);
    %compute:
    prm.fields=load(field_param_file_name);
    FE_PSTH = FE_compute_PSTH(FE,prm);

    % with 2D bins size:    
    [PSTH_for_2D_bin,spike_density,time_spent] = computePSTH([bsp_x_pos{:}],frames_per_second,[spikes_x_pos{:}],allo_X_bins_vector_2D,solo_time_spent_minimum_for_1D_bins,ker_SD);

    %% Field detection
    if ~isempty(not_nan_spikes_x_pos)
        %create data for field detection:
        FR_map.PSTH=FE_PSTH.PSTH;
        FR_map.bin_centers=solo_X_bins_vector_of_centers;
        FR_map.bin_size=prm.fields.bin_size;
        %field detection code:
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
    solo(ii_dir).FE_struct_for_tamirs_code=FE;
    solo(ii_dir).PSTH_for_2D_bin=PSTH_for_2D_bin;
    solo(ii_dir).mean_solo_y_per_x_pos=mean_solo_y_per_x_pos;
    
end

end