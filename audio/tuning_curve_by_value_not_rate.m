function binned_value_smoothed = tuning_curve_by_value_not_rate(spikes_behave_vec,behave_vec,behave_bin_edges,spikes_value_vec,minimum_time_spent,behav_samples_per_sec)
spikes_behave_vec = reshape(spikes_behave_vec,[],1);
spikes_value_vec = reshape(spikes_value_vec,[],1);
grouping_var = discretize(spikes_behave_vec,behave_bin_edges);
[value_in_bins,bins_with_data] = groupsummary(spikes_value_vec,grouping_var,'mean');
binned_value = nan(1,length(behave_bin_edges)-1);
binned_value(bins_with_data) = value_in_bins;
timespent_binned = histcounts(behave_vec, behave_bin_edges) / behav_samples_per_sec;
binned_value(timespent_binned < minimum_time_spent) = NaN;
binned_value_smoothed = fn_nan_smooth_wind_3(binned_value);
binned_value_smoothed = binned_value_smoothed' + (binned_value.*0);
% binned_value_smoothed = reshape(binned_value_smoothed,1,[]);
end

%zeros(1,length(behave_bin_edges)-1)