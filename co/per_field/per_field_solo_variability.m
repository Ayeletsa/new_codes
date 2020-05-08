function solo_firing_rates_time_bins=per_field_solo_variability(field_edges,solo_data,time_X_bin_size,frames_per_second)
% empty array:initiate
solo_firing_rates_time_bins=[];
%find all bsp ind within the field
[bsp_flight_ind,bsp_ind_within_flight]=find(solo_data.bsp.x_pos>field_edges(1) & solo_data.bsp.x_pos<field_edges(2));
flight_vec=unique(bsp_flight_ind);
%run over flights within field and find firing rate within
%a small time windows (same as the bins used in the per
%field by time)
for flight_ii=1:length(flight_vec)
    flight_i=flight_vec(flight_ii);
    relevant_ind=sub2ind(size(solo_data.bsp.x_pos),bsp_flight_ind(find(bsp_flight_ind==flight_i)),bsp_ind_within_flight(find(bsp_flight_ind==flight_i)));
    %find bsp and spike ts withing the field within the
    %flights:
    bsp_relevant_ts=solo_data.bsp.ts_usec(relevant_ind);
    spikes_relevant_ts=solo_data.spikes.ts_usec(find(solo_data.spikes.ts_usec>bsp_relevant_ts(1) & solo_data.spikes.ts_usec<bsp_relevant_ts(end)));
    %create time vector:
    time_vec=bsp_relevant_ts(1):time_X_bin_size:bsp_relevant_ts(end);
    if length(time_vec)>1
        %compute histograms and divide them to get firing rates:
        hist_bsp=hist(bsp_relevant_ts,time_vec);
        hist_spikes=hist(spikes_relevant_ts,time_vec);
        firing_rate_solo_bins=(hist_spikes./hist_bsp)*frames_per_second;
        %add firing rate to all firing rates within this field:
        solo_firing_rates_time_bins=[solo_firing_rates_time_bins,firing_rate_solo_bins];
    end
end

