function co_shuffling_new(co_shuffle_param_file_name)

load(co_shuffle_param_file_name)

files = dir(cell_co_solo_initial_analysis_struct_folder);
behavior_struct_names = {files.name};
%% for each cell
for ii_cell = 3:length(behavior_struct_names)
    ii_cell
    % load cell's egocentric struct
    struct_name = behavior_struct_names{ii_cell};
    file_name = fullfile(cell_co_solo_initial_analysis_struct_folder,struct_name);
    load(file_name);
    behavior_struct=cell_co_solo_initial_analysis;
    
    % save basic information about the cell
    shuffling_struct(1).info.bat = behavior_struct.exp_data.bat;
    shuffling_struct(1).info.day = behavior_struct.exp_data.day;
    shuffling_struct(1).info.cell_num = behavior_struct.exp_data.cell_num;
    
    % shuffle each fircetion separately
    for ii_dir = 1:2
        
        % a. create matrices of shuffled spikes position and distances
        shuffled_data = shuffle_spikes_within_flights (behavior_struct.co(ii_dir).bsp  , behavior_struct.co(ii_dir).spikes,n_shuffles,ego_shuffle);
        
        % b. calculate relevant parameters for every shuffle
        % b.1 create empty variables
        ego_firing_rate_mat = zeros(n_shuffles,length(dis_X_bins_vector_of_centers));
        allo_firing_rate_mat = zeros(n_shuffles,length(pos_X_bins_vector_of_centers));
        max_fr_vec = zeros(1,n_shuffles);
        min_fr_vec = zeros(1,n_shuffles);
        mean_fr_vec = zeros(1,n_shuffles);
        std_fr_vec = zeros(1,n_shuffles);
        median_fr_vec = zeros(1,n_shuffles);
        norm_max_vec = zeros(1,n_shuffles);
        norm_min_vec = zeros(1,n_shuffles);
        range_vec = zeros(1,n_shuffles);
        norm_range_vec = zeros(1,n_shuffles);
        max_min_contrast_index_vec = zeros(1,n_shuffles);
        max_mean_contrast_index_vec = zeros(1,n_shuffles);
        information_per_spike_allo_vec = zeros(1,n_shuffles);
        min_mean_contrast_index_vec = zeros(1,n_shuffles);
        max_med_contrast_index_vec = zeros(1,n_shuffles);
        min_med_contrast_index_vec = zeros(1,n_shuffles);
        range_mean_contrast_index_vec = zeros(1,n_shuffles);
        range_med_contrast_index_vec = zeros(1,n_shuffles);
        information_per_spike_ego_vec = zeros(1,n_shuffles);
        cv_ego_vec = zeros(1,n_shuffles);
        
        % b.2 calculate for evert shuffle
        for ii_shuffle = 1:n_shuffles
            
            % b.2.a. egocentric firing rate
            spikes_dis_m_vec = shuffled_data.spikes_dis_m(ii_shuffle,:);
            bsp_dis_m_vec = behavior_struct.co(ii_dir).bsp.dis_m(:);
            bsp_dis_m_vec = bsp_dis_m_vec(~isnan(bsp_dis_m_vec));
            
             [ego_firing_rate_mat(ii_shuffle,:),spike_density,time_spent] = computePSTH(bsp_dis_m_vec,frames_per_second,spikes_dis_m_vec,dis_X_bins_vector,time_spent_minimum_for_1D_bins,ker_SD);
            [information_per_spike_ego_vec(ii_shuffle), ~] = computeSI(ego_firing_rate_mat(ii_shuffle,:),time_spent);
              %  b.2.c. other egocentric parameters
            %%% when changin the code better to do everything here after
            %%% the loop!!!
            max_fr_vec(ii_shuffle) = nanmax(ego_firing_rate_mat(ii_shuffle,:));
            min_fr_vec(ii_shuffle) = nanmin(ego_firing_rate_mat(ii_shuffle,:));
            mean_fr_vec(ii_shuffle) = nanmean(ego_firing_rate_mat(ii_shuffle,:));
            std_fr_vec(ii_shuffle) = nanstd(ego_firing_rate_mat(ii_shuffle,:));
            median_fr_vec(ii_shuffle) = nanmedian(ego_firing_rate_mat(ii_shuffle,:));
            norm_max_vec(ii_shuffle) = max_fr_vec(ii_shuffle)/mean_fr_vec(ii_shuffle);
            norm_min_vec(ii_shuffle) = min_fr_vec(ii_shuffle)/mean_fr_vec(ii_shuffle);
            range_vec(ii_shuffle) = max_fr_vec(ii_shuffle) - min_fr_vec(ii_shuffle);
            norm_range_vec(ii_shuffle) = (max_fr_vec(ii_shuffle) - min_fr_vec(ii_shuffle)) / mean_fr_vec(ii_shuffle);
            max_min_contrast_index_vec(ii_shuffle) = (max_fr_vec(ii_shuffle) - min_fr_vec(ii_shuffle)) / (max_fr_vec(ii_shuffle) + min_fr_vec(ii_shuffle));
            max_mean_contrast_index_vec(ii_shuffle) = (max_fr_vec(ii_shuffle) - mean_fr_vec(ii_shuffle)) / (max_fr_vec(ii_shuffle) + mean_fr_vec(ii_shuffle));
            min_mean_contrast_index_vec(ii_shuffle) = (mean_fr_vec(ii_shuffle) - min_fr_vec(ii_shuffle)) / (mean_fr_vec(ii_shuffle) + min_fr_vec(ii_shuffle));
            range_mean_contrast_index_vec(ii_shuffle) = (range_vec(ii_shuffle) - mean_fr_vec(ii_shuffle)) / (range_vec(ii_shuffle) + mean_fr_vec(ii_shuffle));
            max_med_contrast_index_vec(ii_shuffle) = (max_fr_vec(ii_shuffle) - median_fr_vec(ii_shuffle)) / (max_fr_vec(ii_shuffle) + median_fr_vec(ii_shuffle));
            min_med_contrast_index_vec(ii_shuffle) = (median_fr_vec(ii_shuffle) - min_fr_vec(ii_shuffle)) / (median_fr_vec(ii_shuffle) + min_fr_vec(ii_shuffle));
            range_med_contrast_index_vec(ii_shuffle) = (range_vec(ii_shuffle) - median_fr_vec(ii_shuffle)) / (range_vec(ii_shuffle) + median_fr_vec(ii_shuffle));
            
            %CV
            cv_ego_vec(ii_shuffle)=std_fr_vec(ii_shuffle)/mean_fr_vec(ii_shuffle);
            
            %  b.2.d. allocentric firing rate
            spikes_x_pos_vec = shuffled_data.spikes_x_pos(ii_shuffle,:);
            bsp_x_pos_vec = behavior_struct.co(ii_dir).bsp.x_pos(:);
            bsp_x_pos_vec = bsp_x_pos_vec(~isnan(bsp_x_pos_vec));
            
            bin_size = allo_bin_size;
            bin_limits = allo_bin_limits;
            bin_edges = bin_limits(1):bin_size:bin_limits(end);
            bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))./2;
           
            
            [allo_firing_rate_mat(ii_shuffle,:),spike_density,allo_time_spent_in_bins] = computePSTH(bsp_x_pos_vec,frames_per_second,spikes_x_pos_vec,bin_edges,time_spent_minimum_for_1D_bins,ker_SD);
            [information_per_spike_allo_vec(ii_shuffle), SI_bits_sec] = computeSI(allo_firing_rate_mat(ii_shuffle,:),allo_time_spent_in_bins);

%             [allo_time_spent_in_bins, ~, ~, allo_firing_rate_mat(ii_shuffle,:), ~,~] ...
%                 = fn_compute_generic_1D_tuning_new_smooth ...
%                 (bsp_x_pos_vec,spikes_x_pos_vec, pos_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
%             
            %  b.2.e. allocentric information
%             fr_nan_ind = isnan(allo_firing_rate_mat(ii_shuffle,:));
%             information_per_spike_allo_vec(ii_shuffle) = spikes_information ( allo_firing_rate_mat(ii_shuffle,~fr_nan_ind), allo_time_spent_in_bins(~fr_nan_ind));
%             
        end
        
        %%%%  TO DO!! 
        %%% also I don't like the word params for data, these are not
        %%% parameters, confusing!
        
        % cv calculations
        params.cv.values=cv_ego_vec;
        
        
        % c. save them into 'params' struct
        params.max.values = max_fr_vec;
        params.min.values = min_fr_vec;
        params.mean.values = mean_fr_vec;
        params.std.values = std_fr_vec;
        params.median.values = median_fr_vec;
        params.norm_max.values =norm_max_vec;
        params.norm_min.values =norm_min_vec;
        params.range.values = range_vec;
        params.norm_range.values = norm_range_vec;
        params.max_min_contrast_index.values = max_min_contrast_index_vec;
        params.max_mean_contrast_index.values = max_mean_contrast_index_vec;
        params.min_mean_contrast_index.values = min_mean_contrast_index_vec;
        params.range_mean_contrast_index.values = range_mean_contrast_index_vec;
        params.max_med_contrast_index.values = max_med_contrast_index_vec;
        params.min_med_contrast_index.values = min_med_contrast_index_vec;
        params.range_med_contrast_index.values = range_med_contrast_index_vec;
        params.information_per_spike_ego.values = information_per_spike_ego_vec;
        params.information_per_spike_allo.values = information_per_spike_allo_vec;
        
        % d. calculate significance for each parameter
        params = significance_for_fr_params(params);
        
        
        
        % zscore calculation:
        zscore_shuffle=zscore(ego_firing_rate_mat,0);
        shuffled_data.zscore_of_tuning_curve_data=zscore_shuffle(1,:);
        shuffled_data.zscore_of_shuffle=zscore_shuffle(2:end,:);
        % e. transfer all data into a single directional struct
        shuffled_data.ego_bin_centers = dis_X_bins_vector_of_centers;
        shuffled_data.allo_bin_centers = pos_X_bins_vector_of_centers;
        shuffled_data.ego_firing_rate = ego_firing_rate_mat;
        shuffled_data.allo_firing_rate = allo_firing_rate_mat;
        shuffled_data.params = params;
        
        
        
        % f. transfer directional data into a single struct
        shuffling_struct(ii_dir).shuffled_data = shuffled_data;
        
        % g. calculate odd-even coherence
        n_co = size(behavior_struct.co(ii_dir).bsp.dis_m,1);
        first_half_ind = 1:2:n_co;
        second_half_ind = 2:2:n_co;
        allo_ego = [0 1]; % egocentric coherence
        
        shuffling_struct(ii_dir).odd_even_coherence = coherence_between_flights...
            (behavior_struct.co(ii_dir).bsp,behavior_struct.co(ii_dir).spikes,first_half_ind,second_half_ind,coherence_X_bins_vector_of_centers,allo_ego,time_spent_minimum_for_1D_bins,frames_per_second);
    end
    
    
    %% save shuffled data struct
    
    if ~exist(co_shuffle_folder_name)
        mkdir(co_shuffle_folder_name)
    end
    shuffle_struct_name = ['co_shuffling_struct_b',num2str( shuffling_struct(1).info.bat),'_d',num2str( shuffling_struct(1).info.day),'_c',num2str(shuffling_struct(1).info.cell_num)];
    file_name=fullfile(co_shuffle_folder_name,shuffle_struct_name);
    save(file_name,'shuffling_struct')
    clear shuffling_struct
end

end