function inclusion=check_if_valid(inclusion,exp_data,ii_dir,param_file_name)
load(param_file_name)
conditions = [];
conditions(1) = exp_data.num_co_per_dir(ii_dir) >= min_co_per_dir;
conditions(2) = exp_data.spikes_num_air(ii_dir)   >=min_n_spike_in_air ;
if isempty(exp_data.mean_fr)
    conditions(3) =0;
else
conditions(3) = exp_data.mean_fr < max_for_pyramidal ;
end
inclusion(ii_dir).conditions_valid_cell = conditions;
inclusion(ii_dir).pyr = conditions(3);
inclusion(ii_dir).valid_cell = all(conditions);
