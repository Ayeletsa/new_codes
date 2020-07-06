function  co = spike_struct_co_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_data,co_param_file_name,per_field_params_file_name,field_param_file_name,co_struct_name)
load(co_param_file_name)
us_factor=1e6;
%% for field detection durin co:
prm.fields=load(field_param_file_name);
prm.fields.min_flights_with_spikes_prc=0.1;
prm.fields.min_flights_with_spikes=3;
prm.fields.min_time_spent_per_bin=0.2;
%NEED TO TAKE IT TO PARAMS IF I WILL USE IT
%
%% load spikes:
spikes_ts = cell_struct.spikes.spikes_ts_usec;
%find spike ts that are not during bsp nan:
spike_pos=interp1(bsp_proc_data(tag_i).ts,bsp_proc_data(tag_i).pos(:,1),spikes_ts);
spikes_ts(isnan(spike_pos))=[];
%% loop for the 2 directions

for ii_dir = 1:2
    %load co behavior struct:
    co_dir_bsp_struct_to_load=[co_struct_name,'_dir_',num2str(ii_dir),'.mat'];
    load(co_dir_bsp_struct_to_load)
    %% remove flights for unstabel cells:
    [logical_vec_of_active_flight]=correct_behavior_to_work_only_for_active_time_per_cell(cell_struct,bsp_ts_usec);
    ind_of_relevant_flight=find(logical_vec_of_active_flight);

    bsp_ts_usec=bsp_ts_usec(ind_of_relevant_flight);
    bsp_dis_m_at_co=bsp_dis_m_at_co(ind_of_relevant_flight);
    bsp_time_to_co=bsp_time_to_co(ind_of_relevant_flight);
    bsp_vel_x_co=bsp_vel_x_co(ind_of_relevant_flight);
    bsp_vel_xy_at_co=bsp_vel_xy_at_co(ind_of_relevant_flight);
    bsp_x_pos_at_co=bsp_x_pos_at_co(ind_of_relevant_flight);
    bsp_y_diff_co=bsp_y_diff_co(ind_of_relevant_flight);
    bsp_y_pos_at_co=bsp_y_pos_at_co(ind_of_relevant_flight);
    all_dir_co_times_usec=all_dir_co_times_usec(ind_of_relevant_flight);
    
    bsp.ind_of_relevant_flight = ind_of_relevant_flight;

%%

    n_co_points=length(all_dir_co_times_usec);

    % initiate spikes cells:
    spikes_ts_usec = cell(n_co_points,1);
    spikes_time_to_co = cell(n_co_points,1);
    spikes_dis_m_at_co = cell(n_co_points,1);
    spikes_x_pos_at_co = cell(n_co_points,1);
    spikes_y_pos_at_co = cell(n_co_points,1);
    spikes_y_diff_co = cell(n_co_points,1);
    
    % loop for every cross over
    for ii_co = 1:n_co_points
        
         % find co time
         co_time_usec = all_dir_co_times_usec(ii_co);

        % a. find spikes index
        spikes_relative_times = spikes_ts - co_time_usec;      
        spikes_ind = find(spikes_ts> bsp_ts_usec{ii_co}(1) & spikes_ts<bsp_ts_usec{ii_co}(end));
   
        % b. find spikes values
        spikes_ts_usec{ii_co} = spikes_ts(spikes_ind);
        spikes_time_to_co{ii_co} = spikes_relative_times(spikes_ind);
        
        spikes_dis_m_at_co{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_dis_m_at_co{ii_co},spikes_ts_usec{ii_co});
        spikes_x_pos_at_co{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_x_pos_at_co{ii_co},spikes_ts_usec{ii_co});
        spikes_y_pos_at_co{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_y_pos_at_co{ii_co},spikes_ts_usec{ii_co});
        spikes_y_diff_co{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_y_diff_co{ii_co},spikes_ts_usec{ii_co});
        
%         % if direction is negative, flip relative positions so negative locations will represent positions before the CO
%         if ii_dir == 2
%             spikes_dis_m_at_co{ii_co} = - spikes_dis_m_at_co{ii_co};
%             spikes_y_diff_co{ii_co} = -  spikes_y_diff_co{ii_co};
%         end
%         
    end
    
    
    %% turn cells into mat padded with NaNs
    
    % bsp
    maxSize = max(cellfun(@numel,bsp_time_to_co));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    % a. time to co
    rmat = cellfun(fcn,bsp_time_to_co,'UniformOutput',false);
    bsp.time_to_co = vertcat(rmat{:});
    % b. relative distance
    rmat = cellfun(fcn,bsp_dis_m_at_co,'UniformOutput',false);
    bsp.dis_m = vertcat(rmat{:});
    %c. x position
    rmat = cellfun(fcn,bsp_x_pos_at_co,'UniformOutput',false);
    bsp.x_pos = vertcat(rmat{:});
    %d. y position
    rmat = cellfun(fcn,bsp_y_pos_at_co,'UniformOutput',false);
    bsp.y_pos = vertcat(rmat{:});
    %e. y diff between bats
    rmat = cellfun(fcn,bsp_y_diff_co,'UniformOutput',false);
    bsp.y_diff = vertcat(rmat{:});
    %f. absolute time stamps
    rmat = cellfun(fcn,bsp_ts_usec,'UniformOutput',false);
    bsp.ts_usec =  vertcat(rmat{:});
    
    % spikes
    maxSize = max(cellfun(@numel,spikes_time_to_co));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    % a. time to co
    % spikes_time_from_co(cellfun(@isempty,spikes_time_from_co)) = {nan};
    rmat = cellfun(fcn,spikes_time_to_co,'UniformOutput',false);
    spikes.time_to_co = vertcat(rmat{:});
    % b. relative distance
    % spikes_dis_m_at_co(cellfun(@isempty,spikes_dis_m_at_co)) = {nan};
    rmat = cellfun(fcn,spikes_dis_m_at_co,'UniformOutput',false);
    spikes.dis_m = vertcat(rmat{:});
    %c. x position
    % spikes_x_pos_at_co(cellfun(@isempty,spikes_x_pos_at_co)) = {nan};
    rmat = cellfun(fcn,spikes_x_pos_at_co,'UniformOutput',false);
    spikes.x_pos = vertcat(rmat{:});
    %d. y position
    % spikes_y_pos_at_co(cellfun(@isempty,spikes_y_pos_at_co)) = {nan};
    rmat = cellfun(fcn,spikes_y_pos_at_co,'UniformOutput',false);
    spikes.y_pos = vertcat(rmat{:});
    %e. y diff between bats
    % spikes_y_diff_co(cellfun(@isempty,spikes_y_diff_co)) = {nan};
    rmat = cellfun(fcn,spikes_y_diff_co,'UniformOutput',false);
    spikes.y_diff = vertcat(rmat{:});
    %f. absolute time stamps
    % spikes_ts_usec(cellfun(@isempty,spikes_y_diff_co)) = {nan};
    rmat = cellfun(fcn,spikes_ts_usec,'UniformOutput',false);
    spikes.ts_usec =  vertcat(rmat{:});
    
    % save into struct
    co(ii_dir).spikes = spikes;
    co(ii_dir).bsp = bsp;
    
    
    %% calculate basic paramters:
    
    % a. number of CO in each direction
    info.n_co = n_co_points;
    
    % b. n spikes in each direction
    info.n_spikes = sum(sum(isfinite(spikes.x_pos)));
    
    % calculate velocity reduction in the center compared to median
    info.vel=find_vel_prop_reduction(bsp_dis_m_at_co,bsp_vel_xy_at_co,long_dis_thresh,short_dis_thresh);
    
    co(ii_dir).info = info;
    %% calculate tuning curves
    
    % a. time to co
    bsp_vec = bsp.time_to_co(isfinite(bsp.time_to_co));
    spikes_vec = spikes.time_to_co(isfinite(spikes.time_to_co));
    [time_spent, ~, ~, PSTH, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
    
    %[PSTH,spike_density,time_spent] = computePSTH(bsp_vec,frames_per_second,spikes_vec,time_X_bins_vector,time_spent_minimum_for_1D_bins,ker_SD);

    
    firing_rate.time_to_co = [PSTH;time_X_bins_vector_of_centers];
    
    % b. relative distance
    bsp_vec = bsp.dis_m(isfinite(bsp.dis_m));
    spikes_vec = spikes.dis_m(isfinite(spikes.dis_m));
    [time_spent, ~, ~, PSTH, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
   
    %[PSTH,spike_density,time_spent] = computePSTH(bsp_vec,frames_per_second,spikes_vec,dis_X_bins_vector,time_spent_minimum_for_1D_bins,ker_SD);
    [SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent);
    firing_rate.dis_m = [PSTH;dis_X_bins_vector_of_centers];
    firing_rate.dis_m_SI=SI_bits_spike;
    % c. allocentric position
   
    % calculate solo firing rate without overlapping spikes
    data=[bsp_x_pos_at_co';bsp_ts_usec';spikes_x_pos_at_co';spikes_ts_usec'];
    field_names={'pos','ts','spikes_pos','spikes_ts'};
    FE=cell2struct(data,field_names,1);
    FE_PSTH = FE_compute_PSTH(FE,prm);
    allo_co_x_pos_fr=FE_PSTH.PSTH;
    
    firing_rate.allo_x_pos = [allo_co_x_pos_fr;prm.fields.bin_centers];
    % compute field detection on allo during co:
    FR_map.PSTH=FE_PSTH.PSTH;
    FR_map.bin_centers=prm.fields.bin_centers;
    FR_map.bin_size=prm.fields.bin_size;
    
    fields=detect_field_based_on_tamir(FR_map,FE,prm) ;
    co(ii_dir).allo_fields_during_co=fields;
    %
    %
    
    % d. compute 2D tuning allocentric*CO distance between bats
    
    
    %create bin vector:
    
    %data:
    bsp_vec_allo = bsp.x_pos(isfinite(bsp.x_pos));
    spikes_vec_allo = spikes.x_pos(isfinite(spikes.x_pos));
    bsp_vec_dis = bsp.dis_m(isfinite(bsp.dis_m));
    spikes_vec_dis = spikes.dis_m(isfinite(spikes.dis_m));
    if ~isempty(spikes_vec_allo)
    [~, ~, firing_rate.field_density_smoothed_XY_with_NaN, ~, ~, ~, ~, ~, ~, ~] ...
        = fn_compute_generic_2D_field ...
        (dis_X_bins_vector_2D, allo_X_bins_vector_2D, bsp_vec_dis, bsp_vec_allo, spikes_vec_dis, spikes_vec_allo, time_spent_minimum_for_2D_bins, frames_per_second, sigma_a, hsize, legalize_by_neighbor_bins_flag);
    else
        firing_rate.field_density_smoothed_XY_with_NaN=nan*zeros(length(dis_X_bins_vector_2D)-1,length(allo_X_bins_vector_2D)-1);
    end
    
    %allo for the 2D plot -
    bsp_vec = bsp.x_pos(isfinite(bsp.x_pos));
    spikes_vec = spikes.x_pos(isfinite(spikes.x_pos));
    
    [PSTH,spike_density,time_spent] = computePSTH(bsp_vec,frames_per_second,spikes_vec,allo_X_bins_vector_2D,time_spent_minimum_for_1D_bins,ker_SD);

    firing_rate.allo_x_pos_fr_for_2D = PSTH;
    
    %ego (dis) for the 2D plot
    bsp_vec = bsp.dis_m(isfinite(bsp.dis_m));
    spikes_vec = spikes.dis_m(isfinite(spikes.dis_m));
       [~, ~, ~, PSTH, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,dis_X_bins_vector_of_centers_2D, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
   
    
    firing_rate.dis_x_pos_fr_for_2D = PSTH;
    
    
    
    %% per field analysis - calculate egocentric tuning within a place filed:
    fields=solo_data(ii_dir).fields;
    if ~isempty(fields)
        %1. run per field analysis on width based on half hight:
        %======================================================
        field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2);
        [per_field_href,per_field_tunings_corrs]=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
        
%         %2. run per field analysis on width based on prctl:
%         %======================================================
%         field_edges=reshape([fields.edges_prc],2,length([fields.edges_prc])/2);
%         per_field_prc=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
%         
%         %3. run per field analysis on width based on all spikes:
%         %======================================================
%         field_edges=reshape([fields.edges_all_spk],2,length([fields.edges_all_spk])/2);
%         per_field_all_spk=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
%         
        
    else
        per_field_href=struct();
        per_field_tunings_corrs=[];
%         per_field_prc=struct();
%         per_field_all_spk=struct();
    end
    
    
    %save all firing rates data into main struct
    co(ii_dir).firing_rate = firing_rate;
    co(ii_dir).per_field_href=per_field_href;
    co(ii_dir).per_field_tunings_corrs=per_field_tunings_corrs;
%     co(ii_dir).per_field_prc=per_field_prc;
%     co(ii_dir).per_field_all_spk=per_field_all_spk;
    %% Inter field analysis:
    %1. Find if there are fields in co and not in solo:
    %-------------------------------------------------------
    if ~isempty(co(ii_dir).allo_fields_during_co)
        peak_allo_co_loc=[co(ii_dir).allo_fields_during_co.loc];
        solo_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2);
        
        for f_i=1:length(peak_allo_co_loc)
            de_novo_fields=isempty(find(peak_allo_co_loc(f_i)>solo_edges(1,:) & peak_allo_co_loc(f_i)<solo_edges(2,:)));
            co(ii_dir).allo_fields_during_co(f_i).de_novo_fields=de_novo_fields;
        end
        
    end
    %2. find inter fields:
    %-------------------------------------------------------
    
    if ~isempty(fields)
        field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2); %for now run on href
        new_field_edges=[];
        inter_field_edge=[];
        %a. enlarge field size by 50% to each side
        field_sizes=[fields.width_href];
        new_field_edges(1,:)=field_edges(1,:)-0.5.*field_sizes;
        new_field_edges(2,:)=field_edges(2,:)+0.5.*field_sizes;
        %b. define interfield edges
        inter_field_edge(1,:)=[0, new_field_edges(2,:)];
        inter_field_edge(2,:)=[new_field_edges(1,:), tunnel_end];
        %c. find if fields were merged
        if size(inter_field_edge,2)>1
            inter_field_edge=inter_field_edge(:,inter_field_edge(2,:)>inter_field_edge(1,:));
        end
        %d. run over inter field pos
        if ~isempty(inter_field_edge)
            inter_field=per_field_analysis(inter_field_edge,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
        else
            inter_field=struct();
            
        end
        
        for inter_field_i=1:size(inter_field_edge,2)
             %e. compare solo and co allo tuning within inter field:

            x_pos_bin_limits=[inter_field_edge(1,inter_field_i),inter_field_edge(2,inter_field_i)];
            solo_bsp_x_pos = solo_data(ii_dir).bsp.x_pos;
            solo_spikes_x_pos = solo_data(ii_dir).spikes.x_pos;
            co_bsp_x_pos=bsp.x_pos;
            co_spikes_x_pos=spikes.x_pos;
            solo_co_comparison=allo_co_solo_comparison(x_pos_bin_limits,sig_bins_width*2,solo_bsp_x_pos,solo_spikes_x_pos,co_bsp_x_pos,co_spikes_x_pos,frames_per_second,min_flights_per_bin);
            inter_field(inter_field_i).solo_co_comparison=solo_co_comparison;
            % compute SI and firing rate per inter field
             pos = solo_data(ii_dir).bsp.x_pos(solo_data(ii_dir).bsp.x_pos>x_pos_bin_limits(1) & solo_data(ii_dir).bsp.x_pos<x_pos_bin_limits(2)) ;
             spikes_pos=solo_data(ii_dir).spikes.x_pos(solo_data(ii_dir).spikes.x_pos>x_pos_bin_limits(1) & solo_data(ii_dir).spikes.x_pos<x_pos_bin_limits(2)) ;
             pos_fs=frames_per_second;
             bin_edges=prm.fields.bin_edges;
             min_time_spent_per_bin=prm.fields.min_time_spent_per_bin;
             ker_SD=prm.fields.ker_SD;
            [PSTH,spike_density,time_spent] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent_per_bin,ker_SD);
            [SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent);
            inter_field(inter_field_i).SI_during_solo=SI_bits_spike;
            inter_field(inter_field_i).tuning_during_solo=PSTH;
            inter_field(inter_field_i).solo_fr=pos_fs*(length(spikes_pos)./length(pos));
            inter_field(inter_field_i).inter_field_edge=inter_field_edge;
        end
    else
        inter_field=struct();
    end
    co(ii_dir).inter_field=inter_field;
    
    %%  calculate significante diff between co and solo allocentric representation
    % as for the obstacle - divide the tunnel into bins and compare
    % population of firing rats between CO and Solo
    non_nan_ind=find(isfinite(co(ii_dir).firing_rate.allo_x_pos(1,:))); %this define the area of the analysis
    x_pos_bin_limits=[min(prm.fields.bin_centers(non_nan_ind)), max(prm.fields.bin_centers(non_nan_ind))];
    solo_bsp_x_pos = solo_data(ii_dir).bsp.x_pos;
    solo_spikes_x_pos = solo_data(ii_dir).spikes.x_pos;
    co_bsp_x_pos=bsp.x_pos;
    co_spikes_x_pos=spikes.x_pos;
    solo_co_comparison=allo_co_solo_comparison(x_pos_bin_limits,sig_bins_width,solo_bsp_x_pos,solo_spikes_x_pos,co_bsp_x_pos,co_spikes_x_pos,frames_per_second,min_flights_per_bin);
    co(ii_dir).solo_co_comparison=solo_co_comparison;
    
    
    
    clear firing_rate
end
