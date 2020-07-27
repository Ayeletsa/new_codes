function    ego_signif_cells=signif_ego_cell(shuffle_file_name,cell_co_solo_initial_analysis,dir_i,param_file_name)
load(shuffle_file_name)
load(param_file_name)

n_spikes = cell_co_solo_initial_analysis.co(dir_i).info.n_spikes;
even_odd_coherence = shuffling_struct(dir_i).odd_even_coherence.corr;
ego_inf = shuffling_struct(dir_i).shuffled_data.params.information_per_spike_ego.values(1);
ego_inf_p = shuffling_struct(dir_i).shuffled_data.params.information_per_spike_ego.p_value;
cv_p=shuffling_struct(dir_i).shuffled_data.params.cv.p_value;

% conditions :
a = n_spikes > min_spikes;
b = even_odd_coherence > min_even_odd;
c = ego_inf > min_ego_inf;
d = ego_inf_p < alpha_val;
%d=cv_p<alpha;
e=cell_co_solo_initial_analysis.exp_data.mean_fr<max_for_pyramidal;

if min([a,b,c,d,e]) == 1
    ego_signif_cells=1;
else
    ego_signif_cells=0;
end
