function  co = spike_struct_co_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_data,co_param_file_name,per_field_params_file_name,field_param_file_name,co_struct_name)
load(co_param_file_name)
us_factor=1e6;
%% for field detection durin co:
%NEED TO TAKE IT TO PARAMS IF I WILL USE IT

prm.fields=load(field_param_file_name);
prm.fields.min_flights_with_spikes_prc=0.1;
prm.fields.min_flights_with_spikes=3;
prm.fields.min_time_spent_per_bin=0.2;

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
    co(ii_dir).spikes_full = spikes;
    co(ii_dir).bsp_full = bsp;
    
    
    %% calculate tuning curves
    allo_data=[bsp_x_pos_at_co';bsp_ts_usec';spikes_x_pos_at_co';spikes_ts_usec'];
    
    [firing_rate_full,fields]=compute_co_tunings(bsp,spikes,co_param_file_name,allo_data,ii_dir,prm);
    
    co(ii_dir).allo_fields_during_co=fields;
     %% per field analysis - calculate egocentric tuning within a place filed:
        fields=solo_data(ii_dir).fields  ;
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
    
    %% save all firing rates data into main struct
    co(ii_dir).firing_rate_full=firing_rate_full;
    co(ii_dir).per_field_href=per_field_href;
    co(ii_dir).per_field_tunings_corrs=per_field_tunings_corrs;
    co(ii_dir).inter_field=inter_field;
    
    %%  calculate significante diff between co and solo allocentric representation
    % as for the obstacle - divide the tunnel into bins and compare
    % population of firing rats between CO and Solo
    non_nan_ind=find(isfinite(co(ii_dir).firing_rate_full.allo_x_pos(1,:))); %this define the area of the analysis
    x_pos_bin_limits=[min(prm.fields.bin_centers(non_nan_ind)), max(prm.fields.bin_centers(non_nan_ind))];
    solo_bsp_x_pos = solo_data(ii_dir).bsp.x_pos;
    solo_spikes_x_pos = solo_data(ii_dir).spikes.x_pos;
    co_bsp_x_pos=bsp.x_pos;
    co_spikes_x_pos=spikes.x_pos;
    solo_co_comparison=allo_co_solo_comparison(x_pos_bin_limits,sig_bins_width,solo_bsp_x_pos,solo_spikes_x_pos,co_bsp_x_pos,co_spikes_x_pos,frames_per_second,min_flights_per_bin);
    co(ii_dir).solo_co_comparison=solo_co_comparison;
    
    %% compute ego tuning only in rectangle of full coverage:
    
    [bsp,spikes,x_pos_valid]=find_co_data_in_rectangle(firing_rate_full,bsp,spikes,allo_X_bins_vector_2D,allo_X_bins_vector_of_centers_2D,min_coverage_of_ego_bin_2D,combine_rectangle_min_sep);
    
    firing_rate=compute_co_tunings(bsp,spikes,co_param_file_name,allo_data,ii_dir,prm);
    
    
    % save into struct
    co(ii_dir).spikes = spikes;
    co(ii_dir).bsp = bsp;
    co(ii_dir).x_pos_valid_rectangle=x_pos_valid;
    co(ii_dir).firing_rate = firing_rate;

      %% compute corr of ego by allo positions
    if ~isempty(fields)
    field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2);
    width_href=[fields.width_href];
    else
        field_edges=[];
        width_href=[];
    end
    
    if sum(sum(isfinite(spikes.x_pos)))~=0
    count_XY_video=firing_rate.count_XY_video;
    count_XY_spikes=firing_rate.count_XY_spikes;
    
    [corr_ego_by_allo_bin,ind_bin_within_field,allo_distance]=compute_ego_corr_per_allo_bin(width_href,field_edges,count_XY_spikes,count_XY_video,allo_X_bins_vector_of_centers_2D,min_n_spike_ego_corr,dis_X_bins_vector_of_centers_2D,time_spent_minimum_for_2D_bins,frames_per_second);
    co(ii_dir).corr_ego_by_allo_bin.corr_ego_by_allo_bin=corr_ego_by_allo_bin;
    co(ii_dir).corr_ego_by_allo_bin.ind_bin_within_field=ind_bin_within_field;
    co(ii_dir).corr_ego_by_allo_bin.allo_distance=allo_distance;
    else
        co(ii_dir).corr_ego_by_allo_bin.corr_ego_by_allo_bin= nan;
        co(ii_dir).corr_ego_by_allo_bin.ind_bin_within_field=nan;
        co(ii_dir).corr_ego_by_allo_bin.allo_distance=nan;
    end
    %% calculate basic paramters:
    
    % a. number of CO in each direction
    info.n_co = n_co_points;   
    % b. n spikes in each direction
    info.n_spikes = sum(sum(isfinite(spikes.x_pos)));
  
    % calculate velocity reduction in the center compared to median
    info.vel=behavioral_modes.co_data.co_bsp_data(ii_dir).vel  ;
    
    % calculate y deviation in the center compared to median
    info.y_deviation=behavioral_modes.co_data.co_bsp_data(ii_dir).y_deviation ;
    
    co(ii_dir).info = info;
    %%
    
    clear firing_rate
end
