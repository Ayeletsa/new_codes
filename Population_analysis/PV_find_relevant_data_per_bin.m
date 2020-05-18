function  [FE,data_for_shuffle]=PV_find_relevant_data_per_bin(bsp_data,spike_data,ego_bin_start,ego_bin_end,bin_i,n_co)

%1. find relevant ind in ego bin:
relevant_ind_bsp=(bsp_data.flight_dis>=ego_bin_start(bin_i) & bsp_data.flight_dis<ego_bin_end(bin_i));
relevant_ind_spike=(spike_data.flight_dis>=ego_bin_start(bin_i) & spike_data.flight_dis<ego_bin_end(bin_i));

%2. initialize nan mats:
co_i_bsp_pos=nan*zeros(size(bsp_data.flight_x_pos));
co_i_spike_pos=nan*zeros(size(spike_data.flight_x_pos));
co_i_bsp_ts=nan*zeros(size(bsp_data.flight_ts));
co_i_spike_ts=nan*zeros(size(spike_data.flight_ts));

%3. take only relevant data for bin:
co_i_bsp_pos(relevant_ind_bsp)=bsp_data.flight_x_pos(relevant_ind_bsp);
co_i_spike_pos(relevant_ind_spike)=spike_data.flight_x_pos(relevant_ind_spike);
co_i_bsp_ts(relevant_ind_bsp)=bsp_data.flight_ts(relevant_ind_bsp);
co_i_spike_ts(relevant_ind_spike)=spike_data.flight_ts(relevant_ind_spike);
% 
%4. convert to cell arrays:
bsp_x_pos=mat2cell(co_i_bsp_pos,ones(size(co_i_bsp_pos,1),1),[size(co_i_bsp_pos,2)]);
bsp_ts_usec=mat2cell(co_i_bsp_ts,ones(size(co_i_bsp_ts,1),1),[size(co_i_bsp_ts,2)]);
spikes_x_pos=mat2cell(co_i_spike_pos,ones(size(co_i_spike_pos,1),1),[size(co_i_spike_pos,2)]);
spikes_ts_usec=mat2cell(co_i_spike_ts,ones(size(co_i_spike_ts,1),1),[size(co_i_spike_ts,2)]);

%5. create flight struct:
data_for_FE=[bsp_x_pos';bsp_ts_usec';spikes_x_pos';spikes_ts_usec'];
field_names={'pos','ts','spikes_pos','spikes_ts'};
FE=cell2struct(data_for_FE,field_names,1);

% Data needed for shuffle analysis 
data_for_shuffle.relevant_shuffle_size_per_flight=round(sum(relevant_ind_bsp(:))/n_co);
data_for_shuffle.relevant_ind_bsp=relevant_ind_bsp;

