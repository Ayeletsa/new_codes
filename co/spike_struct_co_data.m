function  co = spike_struct_co_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_data,co_param_file_name,per_field_params_file_name,field_param_file_name)
load(co_param_file_name)
us_factor=1e6;
%% for field detection durin co:
prm.fields=load(field_param_file_name);
prm.fields.min_flights_with_spikes_prc=0.1;
prm.fields.min_flights_with_spikes=3;
params.fields.min_time_spent_per_meter=0.2;
%NEED TO TAKE IT TO PARAMS IF I WILL USE IT
%% create vectors of variables


% both bats are flying:
FE_self=bsp_proc_data(tag_i).flight_ind;
FE_other=bsp_proc_data(3-tag_i).flight_ind;
ind_both_bats_flying=intersect(FE_self,FE_other);
% find ind that bats are flying in different directions
temp_bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying,1)';
temp_bsp_other_x_pos = bsp_proc_data(3-tag_i).pos(ind_both_bats_flying,1)';
dir_flight=sign(diff(temp_bsp_x_pos));
dir_flight_other=sign(diff(temp_bsp_other_x_pos));
ind_different_direction=find((dir_flight-dir_flight_other)~=0)+1;
ind_both_bats_flying_different_dirs=ind_both_bats_flying(ind_different_direction);
% a. time stamps
bsp_ts = bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs);
spikes_ts = cell_struct.spikes.spikes_ts_usec;

% b. x position
bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1)';
other_bat_x_pos=bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1);
%spikes_x_pos = interp1(bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs) , bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1),spikes_ts);


% relative distance from the other bat
bsp_dis_m = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1)' - bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1)';
%other_bat_pos_at_spikes = interp1(bsp_proc_data(3-tag_i).ts(ind_both_bats_flying_different_dirs) , bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1),spikes_ts);
%spikes_dis_m = spikes_x_pos - other_bat_pos_at_spikes;

% y position
self_y_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,2)';
%spikes_y_pos = interp1(bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs), self_y_pos,spikes_ts);

% y difference between bats
other_y_pos = bsp_proc_data(3 - tag_i).pos(ind_both_bats_flying_different_dirs,2)';
y_diff = self_y_pos - other_y_pos;
%spikes_y_diff = interp1(bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs), y_diff,spikes_ts);

% cross overs parameters
co_x_positions = bsp_proc_data(tag_i).pos(behavioral_modes.CO_point,1);
co_times_usec = bsp_proc_data(tag_i).ts(behavioral_modes.CO_point);
co_directions = sign(co_x_positions - bsp_proc_data(tag_i).pos(behavioral_modes.CO_point-1,1));
direction_ind = {find(co_directions>0),find(co_directions<0)};


%% loop for the 2 directions

for ii_dir = 1:2
    
    dir_ind = direction_ind{ii_dir};
    n_co_points = length(dir_ind);
    
    % create empty arrays for CO loop
    % a. bsp
    bsp_ts_usec = cell(n_co_points,1);
    bsp_time_to_co = cell(n_co_points,1);
    bsp_dis_m_at_co = cell(n_co_points,1);
    bsp_x_pos_at_co = cell(n_co_points,1);
    bsp_y_pos_at_co = cell(n_co_points,1);
    bsp_y_diff_co = cell(n_co_points,1);
    % b. spikes
    spikes_ts_usec = cell(n_co_points,1);
    spikes_time_to_co = cell(n_co_points,1);
    spikes_dis_m_at_co = cell(n_co_points,1);
    spikes_x_pos_at_co = cell(n_co_points,1);
    spikes_y_pos_at_co = cell(n_co_points,1);
    spikes_y_diff_co = cell(n_co_points,1);
    
    % loop for every cross over
    for ii_co = 1:n_co_points
        
        % find co time
        co_time_usec = co_times_usec(dir_ind(ii_co));
        
        % find bsp samples within a smooth_window of -+time_before_after from CO
        % a. find samples index
        bsp_relative_times = bsp_ts - co_time_usec;
        bsp_time_criteria_ind = ( - time_before_after_co*us_factor < bsp_relative_times & bsp_relative_times < time_before_after_co*us_factor);
        bsp_dis_criteria_ind = ( abs(bsp_dis_m)<dis_before_after_co)';
        bsp_ind = bsp_time_criteria_ind & bsp_dis_criteria_ind;
        
        % b. find samples values
        bsp_ts_usec{ii_co} =  bsp_ts(bsp_ind)';
        bsp_time_to_co{ii_co} = bsp_relative_times(bsp_ind)';
        bsp_dis_m_at_co{ii_co} = bsp_dis_m(bsp_ind);
        bsp_x_pos_at_co{ii_co} = bsp_x_pos(bsp_ind);
        bsp_y_pos_at_co{ii_co} =self_y_pos(bsp_ind);
        bsp_y_diff_co{ii_co} = y_diff(bsp_ind);
        
        % find spikes within a smooth_window of -+time_before_after from CO
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
        
        % if direction is negative, flip relative positions so negative locations will represent positions before the CO
        if ii_dir == 2
            spikes_dis_m_at_co{ii_co} = - spikes_dis_m_at_co{ii_co};
            bsp_dis_m_at_co{ii_co} = - bsp_dis_m_at_co{ii_co};
            spikes_y_diff_co{ii_co} = -  spikes_y_diff_co{ii_co};
            bsp_y_diff_co{ii_co} = - bsp_y_diff_co{ii_co};
        end
        
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
    
    % c. time spent in wide bins
    % the 2 dimentional ego*allo plot is divided into 8 by 8 matrix, and time spent in each bin is measured.
    % This is in order to quantify the 2d coverage
    
    
    time_spent_ind = zeros(1,n_time_spent_bins^2);
    
    for ii_ego = 1:n_time_spent_bins
        ego_limits = [ego_bins_edges(ii_ego),ego_bins_edges(ii_ego+1)];
        ego_ind = (ego_limits(1) <= bsp.dis_m & bsp.dis_m < ego_limits(2));
        for ii_allo = 1:n_time_spent_bins
            allo_limits = [allo_bins_edges(ii_allo),allo_bins_edges(ii_allo+1)];
            allo_ind =  (allo_limits(1) > bsp.x_pos & bsp.x_pos > allo_limits(2));
            ind = (ii_ego-1)*n_time_spent_bins + ii_allo;
            time_spent_ind(ind) = sum(sum(ego_ind & allo_ind))/frames_per_second >time_spent_criteria;
        end
    end
    
    info.time_spent_in_bins =  time_spent_ind; % n_time_spent_bins^2 long index, 1 is for sufficiant time spent;
    
    % this is how to draw this data:
    % N =sqrt(length(time_spent_ind));
    % colors = [.7 .7 .7; 1 1 1];
    % color_inds = time_spent_ind+1;
    % r = colors(color_inds,1);
    % g = colors(color_inds,2);
    % b = colors(color_inds,3);
    % checkers = cat(2,r,g,b);
    % checkers = reshape(checkers,[N,N,3]);
    % imagesc(checkers);
    % axis equal tight;
    
    co(ii_dir).info = info;
    
    
    %% calculate tuning curves
    
    % a. time to co
    bsp_vec = bsp.time_to_co(isfinite(bsp.time_to_co));
    spikes_vec = spikes.time_to_co(isfinite(spikes.time_to_co));
    
    [~, ~, ~, time_to_co_fr, ~,~] ...
        = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spikes_vec,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
    
    firing_rate.time_to_co = [time_to_co_fr;time_X_bins_vector_of_centers];
    
    % b. relative distance
    bsp_vec = bsp.dis_m(isfinite(bsp.dis_m));
    spikes_vec = spikes.dis_m(isfinite(spikes.dis_m));
    
    [~, ~, ~, dis_m_fr, ~,~] ...
        = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spikes_vec,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
    
    firing_rate.dis_m = [dis_m_fr;dis_X_bins_vector_of_centers];
    
    % c. allocentric position
    bsp_vec = bsp.x_pos(isfinite(bsp.x_pos));
    spikes_vec = spikes.x_pos(isfinite(spikes.x_pos));
    
    
    % calculate solo firing rate without overlapping spikes
    data=[bsp_x_pos_at_co';bsp_ts_usec';spikes_x_pos_at_co';spikes_ts_usec'];
    field_names={'pos','ts','spikes_pos','spikes_ts'};
    FE=cell2struct(data,field_names,1);
    FE_PSTH = FE_compute_PSTH(FE,prm);
    allo_co_x_pos_fr=FE_PSTH.PSTH;
    
%     [~, ~, ~, allo_co_x_pos_fr, ~,~] ...
%         = fn_compute_generic_1D_tuning_new_smooth ...
%         (bsp_vec,spikes_vec,allo_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
%     
    firing_rate.allo_x_pos = [allo_co_x_pos_fr;allo_X_bins_vector_of_centers];
    
    % to prevent any overlap of spikes between CO and solo,
    % find overlapping spikes and remove them from solo data
    
    solo_bsp_ts = solo_data(ii_dir).bsp.ts_usec;
    solo_bsp_x_pos = solo_data(ii_dir).bsp.x_pos;
    solo_spikes_ts = solo_data(ii_dir).spikes.ts_usec;
    solo_spikes_x_pos = solo_data(ii_dir).spikes.x_pos;
    
    overlapping_bsp_ind = ismember(solo_bsp_ts, bsp.ts_usec);
    solo_bsp_x_pos(overlapping_bsp_ind) = nan;
    overlapping_spikes_ind = ismember(solo_spikes_ts,spikes.ts_usec);
    solo_spikes_x_pos(overlapping_spikes_ind) = nan;
%     solo_bsp_ts(overlapping_bsp_ind)=nan;
%     solo_spikes_ts(overlapping_spikes_ind)=nan;
% 
%     % calculate solo firing rate without overlapping spikes
%     data=[solo_bsp_ts';bsp_ts_usec';solo_spikes_x_pos';solo_spikes_ts'];
%     field_names={'pos','ts','spikes_pos','spikes_ts'};
%     FE=cell2struct(data,field_names,1);
%     FE_PSTH = FE_compute_PSTH(FE,prm);
%     solo_x_pos_fr=FE_PSTH.PSTH;
%     
%     
%     bsp_vec = solo_bsp_x_pos(isfinite(solo_bsp_x_pos));
%     spikes_vec = solo_spikes_x_pos(isfinite(solo_spikes_x_pos));
%     
%     [~, ~, ~, solo_x_pos_fr, ~,~] ...
%         = fn_compute_generic_1D_tuning_new_smooth ...
%         (bsp_vec,spikes_vec,allo_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
%     
%     firing_rate.solo_x_pos = [solo_x_pos_fr;allo_X_bins_vector_of_centers];
    
    % d. compute 2D tuning allocentric*CO distance between bats
    
    
        %create bin vector:
       
        %data:
        bsp_vec_allo = bsp.x_pos(isfinite(bsp.x_pos));
        spikes_vec_allo = spikes.x_pos(isfinite(spikes.x_pos));
        bsp_vec_dis = bsp.dis_m(isfinite(bsp.dis_m));
        spikes_vec_dis = spikes.dis_m(isfinite(spikes.dis_m));
        
        [~, ~, firing_rate.field_density_smoothed_XY_with_NaN, ~, ~, ~, ~, ~, ~, ~] ...
            = fn_compute_generic_2D_field ...
            (dis_X_bins_vector_2D, allo_X_bins_vector_2D, bsp_vec_dis, bsp_vec_allo, spikes_vec_dis, spikes_vec_allo, time_spent_minimum_for_2D_bins, frames_per_second, sigma_a, hsize, legalize_by_neighbor_bins_flag);
        
        
        %allo for the 2D plot
        bsp_vec = bsp.x_pos(isfinite(bsp.x_pos));
        spikes_vec = spikes.x_pos(isfinite(spikes.x_pos));
        
        
        [~, ~, ~, allo_x_pos_fr_for_2D, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,allo_X_bins_vector_of_centers_2D, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
        
        firing_rate.allo_x_pos_fr_for_2D = allo_x_pos_fr_for_2D;
        
        %ego (dis) for the 2D plot
        bsp_vec = bsp.dis_m(isfinite(bsp.dis_m));
        spikes_vec = spikes.dis_m(isfinite(spikes.dis_m));
        
        
        [~, ~, ~, dis_x_pos_fr_for_2D, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,dis_X_bins_vector_of_centers_2D, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
        
        firing_rate.dis_x_pos_fr_for_2D = dis_x_pos_fr_for_2D;
        
        % solo with the same bin size:
        % calculate solo firing rate without overlapping spikes
        bsp_vec = solo_bsp_x_pos(isfinite(solo_bsp_x_pos));
        spikes_vec = solo_spikes_x_pos(isfinite(solo_spikes_x_pos));
        
        [~, ~, ~, solo_x_pos_fr, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,allo_X_bins_vector_of_centers_2D, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
        
        firing_rate.solo_x_pos_same_bin_as_2D = solo_x_pos_fr;
        

    
    
    %% per field analysis - calculate egocentric tuning within a place filed:
    fields=solo_data(ii_dir).fields;
    if ~isempty(fields)
    %1. run per field analysis on width based on half hight:
    %======================================================
    field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2);
    per_field_href=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
    
    %2. run per field analysis on width based on prctl:
    %======================================================
    field_edges=reshape([fields.edges_prc],2,length([fields.edges_prc])/2);
    per_field_prc=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
    
    %3. run per field analysis on width based on all spikes:
    %======================================================
    field_edges=reshape([fields.edges_all_spk],2,length([fields.edges_all_spk])/2);
    per_field_all_spk=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
    

    else
        per_field_href=struct();
        per_field_prc=struct();
        per_field_all_spk=struct();
    end
    
    
%     field_edges=solo_data(ii_dir).field_edges;
%     
%     if ~isempty(field_edges)
%         for field_i=1:size(field_edges,2)
%             spikes_vec=spikes.dis_m(find(spikes.x_pos>field_edges(1,field_i) & spikes.x_pos<field_edges(2,field_i)));
%             bsp_vec=bsp.dis_m(find(bsp.x_pos>field_edges(1,field_i) & bsp.x_pos<field_edges(2,field_i)));
%             if isempty(spikes_vec)
%                 per_field.dis_x_fr_per_field{field_i}=nan*ones(1,length(per_field.dis_X_bins_vector_of_centers));
%                 per_field.information_per_spike_per_field{field_i}=nan;
%                 per_field.signif_field{field_i}=nan;
%                 per_field.number_of_spikes_per_field{field_i}=0;
%                 per_field.sparsity{field_i}=nan;
%                 per_field.cv{field_i}=nan;
%                 per_field.modulation_depth{field_i}=nan;
%                 per_field.time_signif_field_new_smooth{field_i}=nan;
%                 per_field.time_fr_per_field_new_smooth_smooth_window{field_i}=nan;
%                 per_field.all_per_field_solo_per_field{field_i}=nan;
% 
%             else
%                 all_firing_rate_solo_per_field=[];
% 
%                 %per field for distance:
%                 %=============================================
%                 old_smooth=0;
%                 smooth_window_vec=[3,5,7];
%                 smooth_type=1;
%                 smooth_tol=1;
%                 for wind_i=1:length(smooth_window_vec)
%                     smooth_window=smooth_window_vec(wind_i);
%                     [timespent_binned_new_smooth, ~, ~, dis_x_fr_per_field_new_smooth_smooth_window(wind_i,:), ~,~] ...
%                         = fn_compute_generic_1D_tuning_new_smooth_new_again ...
%                         (bsp_vec,spikes_vec,per_field.dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0,old_smooth,smooth_window,smooth_type,smooth_tol);
%                     [signif_field_new_smooth{wind_i}]=shuffling_per_field_analysis_new_smooth(bsp_vec,spikes_vec,per_field.dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth);
%                     
%                     [spike_flight_ind,spike_ind_within_flight]=find(solo_data(ii_dir).spikes.x_pos>field_edges(1,field_i) & solo_data(ii_dir).spikes.x_pos<field_edges(2,field_i));
%                     
%                 end
%                 
%                 
%                 per_field.signif_field_new_smooth=signif_field_new_smooth;
%                 per_field.dis_x_fr_per_field_new_smooth_smooth_window{field_i}=dis_x_fr_per_field_new_smooth_smooth_window;
%                 
%                 
%                 %per field solo variability:
%                 %===================================
%                 %find all bsp ind within the field
%                 [bsp_flight_ind,bsp_ind_within_flight]=find(solo_data(ii_dir).bsp.x_pos>field_edges(1,field_i) & solo_data(ii_dir).bsp.x_pos<field_edges(2,field_i));
%                 flight_vec=unique(bsp_flight_ind);
%                 %run over flights within field and find firing rate within
%                 %a small time windows (same as the bins used in the per
%                 %field by time)
%                 for flight_ii=1:length(flight_vec)
%                     flight_i=flight_vec(flight_ii);
%                     relevant_ind=sub2ind(size(solo_data(ii_dir).bsp.x_pos),bsp_flight_ind(find(bsp_flight_ind==flight_i)),bsp_ind_within_flight(find(bsp_flight_ind==flight_i)));
%                     %find bsp and spike ts withing the field within the
%                     %flights:
%                     bsp_relevant_ts=solo_data(ii_dir).bsp.ts_usec(relevant_ind);
%                     spikes_relevant_ts=solo_data(ii_dir).spikes.ts_usec(find(solo_data(ii_dir).spikes.ts_usec>bsp_relevant_ts(1) & solo_data(ii_dir).spikes.ts_usec<bsp_relevant_ts(end)));
%                     %create time vector:
%                     time_vec=bsp_relevant_ts(1):time_X_bin_size:bsp_relevant_ts(end);
%                     if length(time_vec)>1
%                     %compute histograms and divide them to get firing rates:
%                     hist_bsp=hist(bsp_relevant_ts,time_vec);
%                     hist_spikes=hist(spikes_relevant_ts,time_vec);
%                     firing_rate_solo_bins=(hist_spikes./hist_bsp)*frames_per_second;
%                     %add firing rate to all firing rates within this field:
%                     all_firing_rate_solo_per_field=[all_firing_rate_solo_per_field,firing_rate_solo_bins];
%                     end
%                 end
%                 
%                 % DEBUG:old per field analysis do I need it??
%                 [timespent_binned, ~, ~, dis_x_fr_per_field, ~,~] ...
%                     = fn_compute_generic_1D_tuning_new_smooth ...
%                     (bsp_vec,spikes_vec,per_field.dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0);
%                 [signif_field]=shuffling_per_field_analysis(bsp_vec,spikes_vec,per_field.dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val);
%                 
%                 timespent_binned(find(isnan(dis_x_fr_per_field)))=nan;
%                 [information_per_spike] = fn_compute_spatial_info (timespent_binned,dis_x_fr_per_field);
%                 per_field.dis_x_fr_per_field{field_i}=dis_x_fr_per_field;
%                 per_field.information_per_spike_per_field{field_i}=information_per_spike;
%                 per_field.signif_field{field_i}=signif_field;
%                 per_field.number_of_spikes_per_field{field_i}=length(spikes_vec);
%                 per_field.all_firing_rate_solo_per_field{field_i}=all_firing_rate_solo_per_field;
%                 r=dis_x_fr_per_field;
%                 per_field.sparsity{field_i}=(nanmean(r).^2)/nanmean(r.^2);
%                 std_r=nanstd(r);
%                 mean_r=nanmean(r);
%                 per_field.cv{field_i}=std_r/mean_r;
%                 per_field.modulation_depth{field_i}=(max(r)-min(r))./(max(r));
%                 
%                 %per field for time:
%                 %===================================
%                 spikes_vec=spikes.time_to_co(find(spikes.x_pos>field_edges(1,field_i) & spikes.x_pos<field_edges(2,field_i)));
%                 bsp_vec=bsp.time_to_co(find(bsp.x_pos>field_edges(1,field_i) & bsp.x_pos<field_edges(2,field_i)));
%                 smooth_window_vec=[3,5,7];
%                 smooth_type=1;
%                 smooth_tol=1;
%                 for wind_i=1:length(smooth_window_vec)
%                     smooth_window=smooth_window_vec(wind_i);
%                     [timespent_binned_new_smooth, ~, ~, time_fr_per_field_new_smooth_smooth_window(wind_i,:), ~,~] ...
%                         = fn_compute_generic_1D_tuning_new_smooth_new_again ...
%                         (bsp_vec,spikes_vec,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0,old_smooth,smooth_window,smooth_type,smooth_tol);
%                     [time_signif_field_new_smooth{wind_i}]=shuffling_per_field_analysis_new_smooth(bsp_vec,spikes_vec,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth);
%                     
%                 end
%                 per_field.time_signif_field_new_smooth{field_i}=time_signif_field_new_smooth;
%                 per_field.time_fr_per_field_new_smooth_smooth_window{field_i}=time_fr_per_field_new_smooth_smooth_window;
%                 
%                 %                 [time_timespent_binned, ~, ~, time_fr_per_field, ~,~] ...
%                 %                     = fn_compute_generic_1D_tuning_new_smooth ...
%                 %                     (bsp_vec,spikes_vec,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0);
%                 %                 [time_signif_field]=shuffling_per_field_analysis(bsp_vec,spikes_vec,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val);
%                 %                 [information_per_spike] = fn_compute_spatial_info (time_timespent_binned,time_fr_per_field);
%                 %                 per_field.time_x_fr_per_field{field_i}=time_fr_per_field;
%                 %                 per_field.time_signif_field{field_i}=time_signif_field;
%                 %                 per_field.time_information_per_spike_per_field{field_i}=information_per_spike;
%             end
%         end
%     
%         
%     end
%     
%save all firing rates data into main struct
co(ii_dir).firing_rate = firing_rate;
co(ii_dir).per_field_href=per_field_href;
co(ii_dir).per_field_prc=per_field_prc;
co(ii_dir).per_field_all_spk=per_field_all_spk;
%% Inter field analysis:

if ~isempty(fields)
    field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2); %for now run on href
    new_field_edges=[];
    inter_field_edge=[];
    %1. enlarge field size by 50% to each side
    field_sizes=[fields.width_href];
    new_field_edges(1,:)=field_edges(1,:)-0.5.*field_sizes;
    new_field_edges(2,:)=field_edges(2,:)+0.5.*field_sizes;
    %2. define interfield edges
    inter_field_edge(1,:)=[0, new_field_edges(2,:)];
    inter_field_edge(2,:)=[new_field_edges(1,:), tunnel_end];
    %3. find if fields were merged
    if size(inter_field_edge,2)>1
        inter_field_edge=inter_field_edge(:,inter_field_edge(2,:)>inter_field_edge(1,:));
    end
    %4. run over inter field pos
    if ~isempty(inter_field_edge)
        inter_field=per_field_analysis(inter_field_edge,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
    else
           inter_field=struct();
 
    end
    
else
    inter_field=struct();
end
co(ii_dir).inter_field=inter_field;
%             inter_field_count=0;
%             for inter_field_i=1:size(inter_field_edge,2)
%                 inter_field_count=inter_field_count+1;
%                 spikes_vec=spikes.dis_m(find(spikes.x_pos>inter_field_edge(1,inter_field_i) & spikes.x_pos<inter_field_edge(2,inter_field_i)));
%                 bsp_vec=bsp.dis_m(find(bsp.x_pos>inter_field_edge(1,inter_field_i) & bsp.x_pos<inter_field_edge(2,inter_field_i)));
%                 
%                 
%                 if isempty(spikes_vec)
%                     
%                     inter_field_firing_rate.dis_x_fr_per_field{inter_field_count}=nan*ones(1,length(firing_rate.dis_X_bins_vector_of_centers));
%                     inter_field_firing_rate.information_per_spike_per_field{inter_field_count}=nan;
%                     inter_field_firing_rate.signif_field{inter_field_count}=nan;
%                     inter_field_firing_rate.number_of_spikes_per_field{inter_field_count}=0;
%                     inter_field_firing_rate.sparsity{inter_field_count}=nan;
%                     inter_field_firing_rate.cv{inter_field_count}=nan;
%                     inter_field_firing_rate.modulation_depth{inter_field_count}=nan;
%                 else
%                     %per field for distance:
%                     [timespent_binned, ~, ~, dis_x_fr_per_field, ~,~] ...
%                         = fn_compute_generic_1D_tuning_new_smooth ...
%                         (bsp_vec,spikes_vec,firing_rate.dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0);
%                     [signif_field]=shuffling_per_field_analysis(bsp_vec,spikes_vec,firing_rate.dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val);
%                     
%                     timespent_binned(find(isnan(dis_x_fr_per_field)))=nan;
%                     [information_per_spike] = fn_compute_spatial_info (timespent_binned,dis_x_fr_per_field);
%                     inter_field_firing_rate.dis_x_fr_per_field{inter_field_count}=dis_x_fr_per_field;
%                     inter_field_firing_rate.information_per_spike_per_field{inter_field_count}=information_per_spike;
%                     inter_field_firing_rate.signif_field{inter_field_count}=signif_field;
%                     inter_field_firing_rate.number_of_spikes_per_field{inter_field_count}=length(spikes_vec);
%                     r=dis_x_fr_per_field;
%                     inter_field_firing_rate.sparsity{inter_field_count}=(nanmean(r).^2)/nanmean(r.^2);
%                     std_r=nanstd(r);
%                     mean_r=nanmean(r);
%                     inter_field_firing_rate.cv{inter_field_count}=std_r/mean_r;
%                     inter_field_firing_rate.modulation_depth{inter_field_count}=(max(r)-min(r))./(max(r));
%                 end
%             end
%             %save all firing rates data into main struct
%             inter_field_firing_rate.inter_field_edge=inter_field_edge;
%             co(ii_dir).inter_field_firing_rate = inter_field_firing_rate;
%         else
%             co(ii_dir).inter_field_firing_rate.inter_field_edge =nan;
%         end
%     else
%         co(ii_dir).inter_field_firing_rate.inter_field_edge =nan;
%     end
%     
    %%  calculate significante diff between co and solo allocentric representation
    % as for the obstacle - divide the tunnel into bins and compare
    % population of firing rats between CO and Solo
    
    % create bins only where you have behavioral coverage during CO
    allo_fr_ind = find(isfinite(co(ii_dir).firing_rate.allo_x_pos(1,:))); % find where we could calculate allocantric firing rate
    allo_fr_limits_ind = [allo_fr_ind(1) allo_fr_ind(end)];
    allo_fr_limits = co(ii_dir).firing_rate.allo_x_pos(2,allo_fr_limits_ind);
    allo_bin_width = mean(diff(co(ii_dir).firing_rate.allo_x_pos(2,:))); % use the same bins width as during allocentric tuning curve calculation
    x_pos_bin_limits = [allo_fr_limits(1)-allo_bin_width/2,allo_fr_limits(2)-allo_bin_width/2];
    % because all bins are in a fixed width, calculate where exactly to begin and end the comparison
    x_pos_mod = mod(diff(x_pos_bin_limits),sig_bins_width);
    sig_bins_limits = x_pos_bin_limits(1)-x_pos_mod/2:sig_bins_width:x_pos_bin_limits(2)+x_pos_mod;
    n_bins = length(sig_bins_limits) - 1;
    
    % create empty variables
    p_fr = zeros(1,n_bins);
    diff_mean_fr = zeros(1,n_bins);
    n_solo_flights = size(solo_bsp_x_pos,1);
    
    % for every bin, check if CO firing rate is significantly different
    % from SOLO firing rate
    for ii_bin = 1:n_bins
        
        bin_start = sig_bins_limits(ii_bin);
        bin_end = sig_bins_limits(ii_bin+1);
        solo_fr_at_bin = zeros(1,n_solo_flights);
        co_fr_at_bin = zeros(1,n_co_points);
        
        %for every solo flight, calculate firing rate at bin
        for ii_solo_flight = 1:n_solo_flights
            solo_bsp_at_flight = solo_bsp_x_pos(ii_solo_flight,:);
            n_bsp_at_bin = sum(bin_start <= solo_bsp_at_flight & solo_bsp_at_flight < bin_end);
            solo_spikes_at_flight = solo_spikes_x_pos(ii_solo_flight,:);
            n_spikes_at_bin = sum(bin_start <= solo_spikes_at_flight & solo_spikes_at_flight < bin_end);
            solo_fr_at_bin(ii_solo_flight) = (n_spikes_at_bin/n_bsp_at_bin) * frames_per_second;
        end
        
        %for every cross-over, calculate firing rate at bin
        for ii_co_flight = 1:n_co_points
            co_bsp_at_flight = bsp.x_pos(ii_co_flight,:);
            n_bsp_at_bin = sum(bin_start <= co_bsp_at_flight & co_bsp_at_flight < bin_end);
            co_spikes_at_flight = spikes.x_pos(ii_co_flight,:);
            n_spikes_at_bin = sum(bin_start <= co_spikes_at_flight & co_spikes_at_flight < bin_end);
            co_fr_at_bin(ii_co_flight) = (n_spikes_at_bin/n_bsp_at_bin) * frames_per_second;
        end
        
        solo_fr_at_bin = solo_fr_at_bin(isfinite(solo_fr_at_bin));
        co_fr_at_bin = co_fr_at_bin(isfinite(co_fr_at_bin));
        % calculate p only for bins with more than n flights
        if length(co_fr_at_bin) >= min_flights_per_bin
            [~,p_fr(ii_bin)] = ttest2(solo_fr_at_bin,co_fr_at_bin);
            diff_mean_fr(ii_bin) = mean(solo_fr_at_bin) - mean(co_fr_at_bin);
        else
            p_fr(ii_bin) = nan;
            diff_mean_fr(ii_bin) = nan;
        end
        
        solo_co_comparison.p_diff_firing_rate = p_fr;
        solo_co_comparison.diff_mean_firing_rate = diff_mean_fr;
        sig_bins_centers = sig_bins_limits(1:end-1) + sig_bins_width/2;
        solo_co_comparison.bin_centers = sig_bins_centers;
        solo_co_comparison.n_bins = n_bins;
        
        co(ii_dir).solo_co_comparison = solo_co_comparison;
        
    end
    
    params.dis_before_after_co = dis_before_after_co;
    params.time_before_after_co = time_before_after_co;
    params.co_positions = {co_x_positions(direction_ind{1}),co_x_positions(direction_ind{2})};
    params.co_ts_usec = {co_times_usec(direction_ind{1}),co_times_usec(direction_ind{2})};
    
    co(1).params = params;
    
    clear firing_rate
end
