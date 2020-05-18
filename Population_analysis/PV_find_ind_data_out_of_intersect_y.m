function ind_intersect=PV_find_ind_data_out_of_intersect_y(ind_intersect,data,range_intersect,allo_X_bins_vector,bin_i,dir_i)


%1. for solo bsp data:
solo_per_bin_bsp_ind=(data(dir_i).solo.bsp.flight_x_pos>=allo_X_bins_vector(bin_i)& data(dir_i).solo.bsp.flight_x_pos<allo_X_bins_vector(bin_i+1));
ind_intersect(bin_i).solo_bsp_ind_for_analysis=zeros(size(data(dir_i).solo.bsp.flight_y_pos));
solo_bsp_y_relevant_ind=(data(dir_i).solo.bsp.flight_y_pos>=range_intersect(1) & data(dir_i).solo.bsp.flight_y_pos<=range_intersect(2));
ind_intersect(bin_i).solo_bsp_ind_for_analysis=(solo_per_bin_bsp_ind==1 & solo_bsp_y_relevant_ind==1);

%2. for solo spike data:
ind_intersect(bin_i).solo_spike_ind_for_analysis=zeros(size(data(dir_i).solo.spikes.flight_y_pos));
solo_x_relevant_spike_ind=(data(dir_i).solo.spikes.flight_x_pos>=allo_X_bins_vector(bin_i)& data(dir_i).solo.spikes.flight_x_pos<allo_X_bins_vector(bin_i+1));
solo_spike_y_relevant_ind=(data(dir_i).solo.spikes.flight_y_pos>=range_intersect(1) & data(dir_i).solo.spikes.flight_y_pos<=range_intersect(2));
ind_intersect(bin_i).solo_spike_ind_for_analysis=(solo_x_relevant_spike_ind==1 & solo_spike_y_relevant_ind==1);

%3. for co bsp data:
ind_intersect(bin_i).co_bsp_ind_for_analysis=zeros(size(data(dir_i).co.bsp.flight_y_pos));
co_x_relevant_bsp_ind=(data(dir_i).co.bsp.flight_x_pos>=allo_X_bins_vector(bin_i)& data(dir_i).co.bsp.flight_x_pos<allo_X_bins_vector(bin_i+1));
co_bsp_y_relevant_ind=(data(dir_i).co.bsp.flight_y_pos>=range_intersect(1) & data(dir_i).co.bsp.flight_y_pos<=range_intersect(2));
ind_intersect(bin_i).co_bsp_ind_for_analysis=(co_x_relevant_bsp_ind==1 & co_bsp_y_relevant_ind==1);


%4. for co spike data:
ind_intersect(bin_i).co_spike_ind_for_analysis=zeros(size(data(dir_i).co.spikes.flight_y_pos));
co_x_relevant_spike_ind=(data(dir_i).co.spikes.flight_x_pos>=allo_X_bins_vector(bin_i)& data(dir_i).co.spikes.flight_x_pos<allo_X_bins_vector(bin_i+1));
co_spike_y_relevant_ind=(data(dir_i).co.spikes.flight_y_pos>=range_intersect(1) & data(dir_i).co.spikes.flight_y_pos<=range_intersect(2));
ind_intersect(bin_i).co_spike_ind_for_analysis=(co_x_relevant_spike_ind==1 & co_spike_y_relevant_ind==1);


