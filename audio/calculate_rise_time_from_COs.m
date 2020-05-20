function ego_tuning_curve = calculate_rise_time_from_COs(bsp,shuffled_spikes_flight,co2use,benf_correction,up_factor,x_bins,time_spent_minimum_for_1D,frames_per_second)
shuffled_spikes_dis = cell2mat(shuffled_spikes_flight(co2use));
n_shuffles = size(shuffled_spikes_dis,1);

bsp2use = bsp(co2use,:);
bsp_vec = bsp2use(~isnan(bsp2use(:)));

%shuffle spikes using only part of the data
ego_shuffle = nan(n_shuffles,length(x_bins));
for shuffle_i = 1:n_shuffles
    spike_vec_shuf = shuffled_spikes_dis(shuffle_i,:);
    
    [~, ~, ~, ego_shuffle(shuffle_i,:), ~,~] = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spike_vec_shuf,x_bins, time_spent_minimum_for_1D, frames_per_second, 0,0,0);
end
ego_tuning_curve.ego_shuffle = ego_shuffle;

%tuning curve peak
[~,ego_tuning_curve.peak_ind] = max(ego_shuffle(1,:));

%calculate significance
prct_alpha_max=prctile(ego_shuffle(2:end,:),100-benf_correction);
prct_alpha_min=prctile(ego_shuffle(2:end,:),benf_correction);
mid_prctile_for_width=prctile(ego_shuffle(2:end,:),50);
sig_pos_bins=find(ego_shuffle(1,:)>prct_alpha_max);
sig_neg_bins=find(ego_shuffle(1,:)<prct_alpha_min);
sig_bins = union(sig_pos_bins,sig_neg_bins);

if ~isempty(sig_bins)
    x_bins_up=linspace(min(x_bins),max(x_bins),length(x_bins)*up_factor);
    non_nan_ind=~isnan(ego_shuffle(1,:));
    ego_true_up=interp1(x_bins(non_nan_ind),ego_shuffle(1,non_nan_ind),x_bins_up);
    non_nan_ind=~isnan(mid_prctile_for_width);
    mid_shuf_up=interp1(x_bins(non_nan_ind),mid_prctile_for_width(non_nan_ind),x_bins_up);
    
    ego_tuning_curve.x_bins_up = x_bins_up;
    ego_tuning_curve.ego_true_up = ego_true_up;
    ego_tuning_curve.mid_shuf_up = mid_shuf_up;
    ego_tuning_curve.width_pos = [];
    ego_tuning_curve.width_neg = [];
    
    if ~isempty(sig_pos_bins)
        lower_than_mid = ego_true_up < mid_shuf_up;
        before_first_sig_bin = x_bins_up <= x_bins(sig_pos_bins(1));
        after_first_sig_bin = x_bins_up >= x_bins(sig_pos_bins(1));
        rise_pos_ind = find(before_first_sig_bin & lower_than_mid,1,'last');
        fall_pos_ind = find(after_first_sig_bin & lower_than_mid,1,'first');
        width_pos = [rise_pos_ind,fall_pos_ind];
        if length(width_pos) == 2
            ego_tuning_curve.width_pos = width_pos;
        end
            
    end
    if ~isempty(sig_neg_bins)
        higher_than_mid = ego_true_up > mid_shuf_up;
        before_first_sig_bin = x_bins_up <= x_bins(sig_neg_bins(1));
        after_first_sig_bin = x_bins_up >= x_bins(sig_neg_bins(1));
        rise_neg_ind = find(before_first_sig_bin & higher_than_mid,1,'last');
        fall_neg_ind = find(after_first_sig_bin & higher_than_mid,1,'first');
        width_neg = [rise_neg_ind,fall_neg_ind];
        if length(width_neg) == 2
            ego_tuning_curve.width_neg = width_neg;
        end
    end
end
    ego_tuning_curve.sig_pos_bins = sig_pos_bins;
    ego_tuning_curve.sig_neg_bins = sig_neg_bins;
    ego_tuning_curve.sig_bins = sig_bins;
end
