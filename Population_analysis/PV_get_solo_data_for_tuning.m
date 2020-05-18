
function FE=PV_get_solo_data_for_tuning(bsp_data,spikes_data)

bsp_x_pos=mat2cell(bsp_data.flight_x_pos,ones(size(bsp_data.flight_x_pos,1),1),[size(bsp_data.flight_x_pos,2)]);
bsp_ts_usec=mat2cell(bsp_data.flight_ts,ones(size(bsp_data.flight_ts,1),1),[size(bsp_data.flight_ts,2)]);
spikes_x_pos=mat2cell(spikes_data.flight_x_pos,ones(size(spikes_data.flight_x_pos,1),1),[size(spikes_data.flight_x_pos,2)]);
spikes_ts_usec=mat2cell(spikes_data.flight_ts,ones(size(spikes_data.flight_ts,1),1),[size(spikes_data.flight_ts,2)]);

data_for_FE=[bsp_x_pos';bsp_ts_usec';spikes_x_pos';spikes_ts_usec'];
field_names={'pos','ts','spikes_pos','spikes_ts'};
FE=cell2struct(data_for_FE,field_names,1);