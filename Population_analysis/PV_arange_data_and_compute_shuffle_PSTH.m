function [PSTH,SI_bits_spike]=PV_arange_data_and_compute_shuffle_PSTH(pos,spikes_pos,prm)

bin_size = prm.fields.bin_size;
bin_limits = prm.fields.bin_limits;
bin_edges = bin_limits(1):bin_size:bin_limits(end);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))./2;
%min_time_spent_per = prm.fields.min_time_spent_per_meter;
min_time_spent_per_bin = prm.fields.min_time_spent_per_bin;
ker_SD = prm.fields.ker_SD;
pos_fs=prm.fields.pos_fs;
[PSTH,spike_density,time_spent] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent_per_bin,ker_SD);
[SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent);
