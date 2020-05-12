function clicks_spikes_analysis(dir_param_file_name,audio_param_file_name,co_param_file_name,solo_param_file_name)
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
for ii_day = 1:length(all_days)
    ii_day
    name_prefix = all_days(ii_day).name;
    clicks_day_struct_name = [clicks_day_struct_folder,'clicks_',name_prefix];
%     if exist(clicks_day_struct_name,'file')
%         continue
%     end
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
    
%     create_clicks_day_struct(p,clicks_day_struct_name,audio_filt_name,audio_param_file_name,general_behavior_data_file_name,tag_i,behav_struct_name,co_param_file_name,solo_param_file_name,click_detection_folder);

    % 3. calculate calling-rate + 4. divide COs by calling-rate rise time
    calling_rate = divide_early_late(clicks_day_struct_name,click_rate_folder,name_prefix,co_param_file_name);
    calling_rate_struct_name = strrep(clicks_day_struct_name,'clicks_','calling_rate_');
    save(calling_rate_struct_name,'calling_rate')
end

return

%% for each cell:
load(co_param_file_name)
cell_files = dir(cells_struct_dir);
spike_files = dir(cell_co_solo_initial_analysis_struct_folder);
shuffle_files = dir(co_shuffle_folder_name);
calling_rate_files = dir(fullfile(clicks_day_struct_folder,'calling_rate_*.mat'));

for ii_cell = 3:length(cell_files)
    ii_cell
    cell_name = cell_files(ii_cell).name;
    
    load(fullfile(cell_files(ii_day).folder,cell_name));
    if ~cell_struct.exp_info.audio_self
        continue
    end
    c_day = num2str(cell_struct.exp_info.day);
    c_bat = num2str(cell_struct.exp_info.bat);
    c_num = num2str(cell_struct.cell_info.cell_num);
    
    calling_rate_struct_name = calling_rate_files(cellfun(@(x) contains(x,c_day) && contains(x,c_bat),{calling_rate_files.name}));
    load(fullfile(calling_rate_struct_name.folder,calling_rate_struct_name.name))
    
    cell_clicks_spikes_analysis.exp_data.bat = c_bat;
    cell_clicks_spikes_analysis.exp_data.day = c_day;
    cell_clicks_spikes_analysis.exp_data.cell_num = c_num;
    
    cell_clicks_spikes_analysis.calling_rate = calling_rate;
    
    
    %% analyze spikes+clicks
    load(fullfile(cell_co_solo_initial_analysis_struct_folder,spike_files(ii_cell).name));
    load(fullfile(co_shuffle_folder_name,shuffle_files(ii_cell).name));
    
    for ii_dir = 1:2
        bsp_dis = cell_co_solo_initial_analysis.co(ii_dir).bsp.dis_m;
        spikes_dis = cell_co_solo_initial_analysis.co(ii_dir).spikes.dis_m;
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
        
        % decide on criterion!!!
        % let's use solo for now
        co_early = cell_clicks_spikes_analysis.calling_rate(ii_dir).solo_early;
        %     co_late = cr(ii_dir).solo_late;
        co_bad = cell_clicks_spikes_analysis.calling_rate(ii_dir).high_bl | cell_clicks_spikes_analysis.calling_rate(ii_dir).solo_prob;
        
        co_early = find(co_early(~co_bad));
        %     co_late = co_late(~co_bad);
        bsp_dis = bsp_dis(~co_bad,:);
        shuffled_spikes_dis_flight = shuffled_spikes_dis_flight(~co_bad);
        nflights = length(shuffled_spikes_dis_flight);
        
        n_permutations = 1000;
        up_factor = 10;
        
        benf_correction=alpha_val/length(dis_X_bins_vector_of_centers); %correct for number of bins
        benf_correction=benf_correction/2; %correct for checking both min and max
        
        perm_indices = zeros(floor(nflights/2),n_permutations);
        delta_sig = nan(1,n_permutations);
        delta_rise = nan(1,n_permutations);
        
        %     addAttachedFiles(gcp,co_param_file_name)
        dis_X_bins_vector_of_centers = dis_X_bins_vector_of_centers;
        time_spent_minimum_for_1D_bins = time_spent_minimum_for_1D_bins;
        frames_per_second = frames_per_second;
        %     parfor ii_perm = 2:n_permutations
        %         perm_ind_early = randsample(nflights,floor(nflights/2));
        %         perm_indices(:,ii_perm) = perm_ind_early;
        %         perm_ind_late = setdiff(1:nflights,perm_ind_early);
        %
        %         ego_early = calculate_rise_time_from_COs ...
        %             (bsp_dis,shuffled_spikes_dis_flight,perm_ind_early,benf_correction,up_factor, ...
        %             dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
        %
        %         ego_late = calculate_rise_time_from_COs ...
        %             (bsp_dis,shuffled_spikes_dis_flight,perm_ind_late,benf_correction,up_factor, ...
        %             dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
        %
        %         if ~isempty(ego_early.rise_dist) && ~isempty(ego_late.rise_dist)
        %             delta_rise(ii_perm) = ego_late.rise_dist - ego_early.rise_dist;
        %         end
        %         if ~isempty(ego_early.first_sig_dist) && ~isempty(ego_late.first_sig_dist)
        %             delta_sig(ii_perm) = ego_late.first_sig_dist - ego_early.first_sig_dist;
        %         end
        %     end
        
        perm_indices(:,1) = co_early;
        perm_ind_late = setdiff(1:nflights,co_early);
        
        ego_early = calculate_rise_time_from_COs ...
            (bsp_dis,shuffled_spikes_dis_flight,co_early,benf_correction,up_factor, ...
            dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
        
        ego_late = calculate_rise_time_from_COs ...
            (bsp_dis,shuffled_spikes_dis_flight,perm_ind_late,benf_correction,up_factor, ...
            dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
        
        if ~isempty(ego_early.rise_dist) && ~isempty(ego_late.rise_dist)
            delta_rise(1) = ego_late.rise_dist - ego_early.rise_dist;
        end
        if ~isempty(ego_early.first_sig_dist) && ~isempty(ego_late.first_sig_dist)
            delta_sig(1) = ego_late.first_sig_dist - ego_early.first_sig_dist;
        end
        
        %     figure
        %     histogram(delta_rise)
        %     figure
        %     histogram(delta_sig)
        %     figure
        %     subplot(1,2,1)
        %     hold on
        %     plot(dis_X_bins_vector_of_centers,ego_early.ego_shuffle(2:end,:),'color',[0.9 0.9 0.9],'LineWidth',3);
        %     plot(dis_X_bins_vector_of_centers,ego_early.ego_shuffle(1,:),'color',[1 0 1],'LineWidth',3);
        %     if ~isempty(ego_early.first_sig_bin) && strcmp(ego_early.sign,'pos')
        %         plot(dis_X_bins_vector_of_centers(ego_early.first_sig_bin),ego_early.ego_shuffle(1,ego_early.first_sig_bin),'g*')
        %     elseif ~isempty(ego_early.first_sig_bin) && strcmp(ego_early.sign,'neg')
        %         plot(dis_X_bins_vector_of_centers(ego_early.first_sig_bin),ego_early.ego_shuffle(1,ego_early.first_sig_bin),'r*')
        %     end
        %     title('early')
        %     subplot(1,2,2)
        %     hold on
        %     plot(dis_X_bins_vector_of_centers,ego_late.ego_shuffle(2:end,:),'color',[0.9 0.9 0.9],'LineWidth',3);
        %     plot(dis_X_bins_vector_of_centers,ego_late.ego_shuffle(1,:),'color',[1 0 1],'LineWidth',3);
        %     if ~isempty(ego_late.first_sig_bin) && strcmp(ego_late.sign,'pos')
        %         plot(dis_X_bins_vector_of_centers(ego_late.first_sig_bin),ego_late.ego_shuffle(1,ego_late.first_sig_bin),'g*')
        %     elseif ~isempty(ego_late.first_sig_bin) && strcmp(ego_late.sign,'neg')
        %         plot(dis_X_bins_vector_of_centers(ego_late.first_sig_bin),ego_late.ego_shuffle(1,ego_late.first_sig_bin),'r*')
        %     end
        %     title('late')
        
        cell_clicks_spikes_analysis.early_late(ii_dir).ego_early = ego_early;
        cell_clicks_spikes_analysis.early_late(ii_dir).ego_late = ego_late;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_rise = delta_rise;
        cell_clicks_spikes_analysis.early_late(ii_dir).delta_sig = delta_sig;
        
        %         bsp = structfun(@(x) x(~high_bl & ~prob,:),bsp,'UniformOutput',false);
        %         spikes = structfun(@(x) x(~high_bl & ~prob,:),spikes,'UniformOutput',false);
        
        
        
        %         figure
        %         clf
        %         p=panel();
        %         c=corr(dis_x_pos_fr_for_2D_early',dis_x_pos_fr_for_2D_late');
        %         sgtitle([cell_name ', dir_' num2str(ii_dir) ', corr=' num2str(c)],'interpreter','none')
        %         p.pack({0.2,[]},2);
        %         %     p.select('all')
        %         %     p.identify();
        %         p(1,1).select()
        %         title('early')
        %         plot(firing_rate.dis_X_bins_vector_of_centers{1,1},dis_x_pos_fr_for_2D_early,'m','LineWidth',2)
        %         p(2,1).select()
        %         fn_plot_2D_field(allo_ego_map_early,firing_rate(ii_dir).dis_X_bins_vector{1},firing_rate(ii_dir).dis_X_bins_vector_of_centers{1},firing_rate(ii_dir).allo_X_bins_vector{1},firing_rate(ii_dir).allo_X_bins_vector_of_centers{1},[]);
        %         axis xy
        %         p(1,2).select()
        %         title('late')
        %         plot(firing_rate.dis_X_bins_vector_of_centers{1,1},dis_x_pos_fr_for_2D_late,'m','LineWidth',2)
        %         p(2,2).select()
        %         fn_plot_2D_field(allo_ego_map_late,firing_rate(ii_dir).dis_X_bins_vector{1},firing_rate(ii_dir).dis_X_bins_vector_of_centers{1},firing_rate(ii_dir).allo_X_bins_vector{1},firing_rate(ii_dir).allo_X_bins_vector_of_centers{1},[]);
        %         axis xy
        %         p.margin = [15 15 15 20];
        
        
    end
    %%
    click_spike_analysis_struct_name = fullfile(clicks_spikes_analysis_struct_folder,['cl_sp_' cell_name]);
    save(click_spike_analysis_struct_name,'cell_clicks_spikes_analysis')
end
end


function ego_tuning_curve = calculate_rise_time_from_COs(bsp,shuffled_spikes_flight,co2use,benf_correction,up_factor,x_bins,time_spent_minimum_for_1D,frames_per_second)
shuffled_spikes_dis = cell2mat(shuffled_spikes_flight(co2use));
n_shuffles = size(shuffled_spikes_dis,1);

bsp2use = bsp(co2use,:);
bsp_vec = bsp2use(~isnan(bsp2use(:)));

ego_shuffle = nan(n_shuffles,length(x_bins));

for shuffle_i = 1:n_shuffles
    spike_vec_shuf = shuffled_spikes_dis(shuffle_i,:);
    
    [~, ~, ~, ego_shuffle(shuffle_i,:), ~,~] = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spike_vec_shuf,x_bins, time_spent_minimum_for_1D, frames_per_second, 0,0,0);
end
prct_alpha_max=prctile(ego_shuffle(2:end,:),100-benf_correction);
prct_alpha_min=prctile(ego_shuffle(2:end,:),benf_correction);
mid_prctile_for_width=prctile(ego_shuffle(2:end,:),50);
first_sig_bin=find(ego_shuffle(1,:)>prct_alpha_max | ego_shuffle(1,:)<prct_alpha_min,1);

ego_tuning_curve.ego_shuffle = ego_shuffle;

if ~isempty(first_sig_bin)
    x_bins_up=linspace(min(x_bins),max(x_bins),length(x_bins)*up_factor);
    non_nan_ind=~isnan(ego_shuffle(1,:));
    ego_true_up=interp1(x_bins(non_nan_ind),ego_shuffle(1,non_nan_ind),x_bins_up);
    non_nan_ind=~isnan(mid_prctile_for_width);
    mid_shuf_up=interp1(x_bins(non_nan_ind),mid_prctile_for_width(non_nan_ind),x_bins_up);
    
    if any(ego_shuffle(1,first_sig_bin)>prct_alpha_max)
        bin_type = 'pos';
        lower_than_mid = find(ego_true_up < mid_shuf_up);
        rise_ind = find(lower_than_mid < first_sig_bin,1,'last');
    elseif any(ego_shuffle(1,first_sig_bin)<prct_alpha_min)
        bin_type = 'neg';
        higher_than_mid = find(ego_true_up > mid_shuf_up);
        rise_ind = find(higher_than_mid < first_sig_bin,1,'last');
    end
    
    ego_tuning_curve.first_sig_dist = x_bins(first_sig_bin);
    ego_tuning_curve.sign = bin_type;
    ego_tuning_curve.rise_dist = x_bins_up(rise_ind);
    
else
    ego_tuning_curve.first_sig_dist = [];
    ego_tuning_curve.rise_dist = [];
end
end