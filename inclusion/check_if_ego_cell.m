function inclusion=check_if_ego_cell(inclusion,co,ii_dir,shuffle,param_file_name)
load(param_file_name)
if inclusion(ii_dir).valid_cell==0
    inclusion(ii_dir).ego_cell=0;
else
% has enough spikes during co
enough_spikes_during_co=co(ii_dir).info.n_spikes>min_spikes_during_co;
% SI > thr
SI_ego=shuffle(ii_dir).shuffled_data.params.information_per_spike_ego.values(1);
SI_thr_signif = SI_ego   > SI_thr_co;
% SI > shuffle
SI_shuffle_signif = SI_ego > ...
    prctile([shuffle(ii_dir).shuffled_data.params.information_per_spike_ego.values(2:end)  ],SI_thr_shuffle_co);
%stability: even vs odd:
r=shuffle(ii_dir).odd_even_coherence.corr;
cell_is_stabel=r>=min_even_odd;
% valid cell
valid_cell=inclusion(ii_dir).valid_cell;
%pyramidal:
pyr=inclusion(ii_dir).pyr;

%% apply all conditions
TF = true;
TF = TF & enough_spikes_during_co;
TF = TF & SI_thr_signif;
TF = TF & SI_shuffle_signif;
TF = TF & cell_is_stabel;
TF = TF & pyr; 
TF = TF & valid_cell; 

inclusion(ii_dir).ego_cell=TF;
end