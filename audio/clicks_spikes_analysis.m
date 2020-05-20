function clicks_spikes_analysis(dir_param_file_name,audio_param_file_name,co_param_file_name,solo_param_file_name,early_late_file_name)
% for each day:
% 1. detect clicks
% 2. assign clicks for each co and solo flight
% 3. calculate calling rate by ego-dis
% 4. divide COs by calling-rate rise time
%
% for each cell:
% 1. early vs. late ego tuning curves

load(dir_param_file_name)

%% for each day:
all_days = dir(fullfile(day_struct_dir,'*.mat'));
for ii_day = 8%:length(all_days)
    ii_day
    name_prefix = all_days(ii_day).name;
    clicks_day_struct_name = [clicks_day_struct_folder,'clicks_',name_prefix];
    
    load(fullfile(all_days(ii_day).folder,name_prefix));
    if ~p.Audio_self
        continue
    end
    
    % 1. detect clicks + 2. assign clicks to COs and solo flights
    tag_self = p.bsp_tag_self;
    tags = [p.bsp_data.tag_ID];
    tag_i = find(ismember(tags,tag_self));
    
    behav_struct_name=[behave_day_struct_folder,'behavioral_modes_',name_prefix];
    general_behavior_data_file_name=[behave_day_struct_folder,'general_behave_analysis_',name_prefix];
    general_behavior_data_file_name = strrep(strrep(general_behavior_data_file_name,'_bat_','_bat'),'_day_','_day');
    audio_filt_name = [audio_struct_dir,strrep(all_days(ii_day).name,'.mat','_audio_filt')];
    
    if ~exist(behav_struct_name,'file')
        continue
    end
    if exist(clicks_day_struct_name,'file')
        load(clicks_day_struct_name)
    else
        create_clicks_day_struct(p,clicks_day_struct_name,audio_filt_name,audio_param_file_name,general_behavior_data_file_name,tag_i,behav_struct_name,co_param_file_name,solo_param_file_name,click_detection_folder);
    end
    %     align_stats
    % 3. calculate calling-rate + 4. divide COs by calling-rate rise time
    %     calling_rate = divide_early_late(clicks_day_struct_name,click_rate_folder,name_prefix,co_param_file_name,audio_param_file_name,early_late_file_name);
    %     calling_rate_struct_name = [clicks_day_struct_folder,'calling_rate_',name_prefix];
    %     save(calling_rate_struct_name,'calling_rate')
end

% return

%% for each cell:
load(co_param_file_name)
load(early_late_file_name)
cell_files = dir(cells_struct_dir);
spike_files = dir(cell_co_solo_initial_analysis_struct_folder);
shuffle_files = dir(co_shuffle_folder_name);
calling_rate_files = dir(fullfile(clicks_day_struct_folder,'calling_rate_*.mat'));

for ii_cell = 80%3:length(cell_files)
    ii_cell
    cell_name = cell_files(ii_cell).name;
    
    load(fullfile(cell_files(ii_day).folder,cell_name));
%     if ~cell_struct.exp_info.audio_self
%         continue
%     end
    c_day = num2str(cell_struct.exp_info.day);
    c_bat = num2str(cell_struct.exp_info.bat);
    c_num = num2str(cell_struct.cell_info.cell_num);
    
    calling_rate_struct_name = calling_rate_files(cellfun(@(x) contains(x,c_day) && contains(x,c_bat),{calling_rate_files.name}));
    if isempty(calling_rate_struct_name)
        continue
    else
        load(fullfile(calling_rate_struct_name.folder,calling_rate_struct_name.name))
    end
    
    cell_clicks_spikes_analysis.exp_data.bat = c_bat;
    cell_clicks_spikes_analysis.exp_data.day = c_day;
    cell_clicks_spikes_analysis.exp_data.cell_num = c_num;
    
    cell_clicks_spikes_analysis.calling_rate = calling_rate;
    
    
    %% analyze spikes+clicks
    spike_struct_name = spike_files(cellfun(@(x) contains(x,c_num) && contains(x,c_day) && contains(x,c_bat) ,{spike_files.name})).name;
    shuff_struct_name = shuffle_files(cellfun(@(x) contains(x,c_num) && contains(x,c_day) && contains(x,c_bat) ,{shuffle_files.name})).name;
    load(fullfile(cell_co_solo_initial_analysis_struct_folder,spike_struct_name));
    load(fullfile(co_shuffle_folder_name,shuff_struct_name));
    
    for ii_dir = 1:2
        bsp_dis = cell_co_solo_initial_analysis.co(ii_dir).bsp.dis_m;
        bsp_allo = cell_co_solo_initial_analysis.co(ii_dir).bsp.dis_m;
        spikes_dis = cell_co_solo_initial_analysis.co(ii_dir).spikes.x_pos;
        spikes_allo = cell_co_solo_initial_analysis.co(ii_dir).spikes.x_pos;
        
        n_spikes_per_flight = sum(~isnan(spikes_dis),2);
        nflights = length(n_spikes_per_flight);
        
        shuffled_spikes_dis = shuffling_struct(ii_dir).shuffled_data.spikes_dis_m;
        n_shuffles = size(shuffled_spikes_dis,1);
        shuffled_spikes_dis_flight = cell(1,nflights);
        
        end_indices = cumsum(n_spikes_per_flight);
        start_indices = end_indices - n_spikes_per_flight + 1;
        
        for ii_flight = 1:nflights
            shuffled_spikes_dis_flight{ii_flight} = shuffled_spikes_dis(:,start_indices(ii_flight):end_indices(ii_flight));
        end
        
        % decide on criterion
        % use half-height criterion
        co_early = cell_clicks_spikes_analysis.calling_rate(ii_dir).half_height_early;
        co_late = cell_clicks_spikes_analysis.calling_rate(ii_dir).half_height_late;
        co_bad = cell_clicks_spikes_analysis.calling_rate(ii_dir).high_bl;
        
        % remove high baseline flights from the data
        co_early(co_bad) = [];
        co_late(co_bad) = [];
        shuffled_spikes_dis_flight(co_bad) = [];
        bsp_dis(co_bad,:) = [];
        bsp_allo(co_bad,:) = [];
        spikes_dis(co_bad,:) = [];
        spikes_allo(co_bad,:) = [];
        nflights = length(shuffled_spikes_dis_flight);
        
        benf_correction=alpha_val/length(dis_X_bins_vector_of_centers); %correct for number of bins
        benf_correction=benf_correction/2; %correct for checking both min and max
        
        % pre-allocate vectors of statistics
        perm_indices = zeros(floor(nflights/2),n_permutations);
        delta_pos_sig = nan(1,n_permutations);
        delta_neg_sig = nan(1,n_permutations);
        delta_pos_rise = nan(1,n_permutations);
        delta_neg_rise = nan(1,n_permutations);
        delta_peak = nan(1,n_permutations);
        
        % Monte Carlo permutation test for ego tuning curves
        dis_X_bins_vector_of_centers = dis_X_bins_vector_of_centers;
        time_spent_minimum_for_1D_bins = time_spent_minimum_for_1D_bins;
        frames_per_second = frames_per_second;
        parfor ii_perm = 2:n_permutations+1
            perm_ind_early = randsample(nflights,floor(nflights/2));
            perm_ind_late = setdiff(1:nflights,perm_ind_early);
            perm_indices(:,ii_perm) = perm_ind_early;
            
            ego_early = calculate_rise_time_from_COs ...
                (bsp_dis,shuffled_spikes_dis_flight,perm_ind_early,benf_correction,tuning_up_factor, ...
                dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
            
            ego_late = calculate_rise_time_from_COs ...
                (bsp_dis,shuffled_spikes_dis_flight,perm_ind_late,benf_correction,tuning_up_factor, ...
                dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
            
            if ~isempty(ego_early.width_pos) && ~isempty(ego_late.width_pos)
                delta_pos_rise(ii_perm) = ego_late.width_pos(1) - ego_early.width_pos(1);
            end
            if ~isempty(ego_early.width_neg) && ~isempty(ego_late.width_neg)
                delta_neg_rise(ii_perm) = ego_late.width_neg(1) - ego_early.width_neg(1);
            end
            if ~isempty(ego_early.sig_pos_bins) && ~isempty(ego_late.sig_pos_bins)
                delta_pos_sig(ii_perm) = ego_late.sig_pos_bins(1) - ego_early.sig_pos_bins(1);
            end
            if ~isempty(ego_early.sig_neg_bins) && ~isempty(ego_late.sig_neg_bins)
                delta_neg_sig(ii_perm) = ego_late.sig_neg_bins(1) - ego_early.sig_neg_bins(1);
            end
            delta_peak(ii_perm) = ego_late.peak_ind - ego_early.peak_ind;
        end
        
        perm_indices(:,1) = find(co_early);
        
        ego_early = calculate_rise_time_from_COs ...
            (bsp_dis,shuffled_spikes_dis_flight,co_early,benf_correction,tuning_up_factor, ...
            dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
        
        ego_late = calculate_rise_time_from_COs ...
            (bsp_dis,shuffled_spikes_dis_flight,co_late,benf_correction,tuning_up_factor, ...
            dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
        
        if ~isempty(ego_early.width_pos) && ~isempty(ego_late.width_pos)
            delta_pos_rise(1) = ego_late.width_pos(1) - ego_early.width_pos(1);
        end
        if ~isempty(ego_early.width_neg) && ~isempty(ego_late.width_neg)
            delta_neg_rise(1) = ego_late.width_neg(1) - ego_early.width_neg(1);
        end
        if ~isempty(ego_early.sig_pos_bins) && ~isempty(ego_late.sig_pos_bins)
            delta_pos_sig(1) = ego_late.sig_pos_bins(1) - ego_early.sig_pos_bins(1);
        end
        if ~isempty(ego_early.sig_neg_bins) && ~isempty(ego_late.sig_neg_bins)
            delta_neg_sig(1) = ego_late.sig_neg_bins(1) - ego_early.sig_neg_bins(1);
        end
        delta_peak(1) = ego_late.peak_ind - ego_early.peak_ind;
        
        bin_size = min(diff(dis_X_bins_vector_of_centers));
        bin_size_up = min(diff(ego_early.x_bins_up));
        
        inv_prcnt = @(x) mean(x(1)<=x(2:end));
        
        
        % ego vs. allo tuning curves
        %early
        bsp_vec_allo = bsp.x_pos(isfinite(bsp.x_pos));
        spikes_vec_allo = spikes.x_pos(isfinite(spikes.x_pos));
        bsp_vec_dis = bsp.dis_m(isfinite(bsp.dis_m));
        spikes_vec_dis = spikes.dis_m(isfinite(spikes.dis_m));
        
        [~, ~, ego_allo_2D_early, ~, ~, ~, ~, ~, ~, ~] = fn_compute_generic_2D_field ...
            (firing_rate.dis_X_bins_vector{bin_dis_i,bin_allo_i}, firing_rate.allo_X_bins_vector{bin_dis_i,bin_allo_i}, bsp_vec_dis, bsp_vec_allo, spikes_vec_dis, spikes_vec_allo, time_spent_minimum_for_2D_bins, frames_per_second, sigma_a, hsize, legalize_by_neighbor_bins_flag);
       
      
        
        %%
        cell_clicks_spikes_analysis.early_late(ii_dir).ego_early = ego_early;
        cell_clicks_spikes_analysis.early_late(ii_dir).ego_late = ego_late;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_pos_sig = delta_pos_sig;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_neg_sig = delta_neg_sig;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_pos_rise = delta_pos_rise;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_neg_rise = delta_neg_rise;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_peak = delta_peak;
        cell_clicks_spikes_analysis.early_late(ii_dir).perm_indices = perm_indices;
    end
    %%
    click_spike_analysis_struct_name = fullfile(clicks_spikes_analysis_struct_folder,['cl_sp_' cell_name]);
    save(click_spike_analysis_struct_name,'cell_clicks_spikes_analysis')
end
end