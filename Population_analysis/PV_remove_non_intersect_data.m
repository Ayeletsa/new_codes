function data=PV_remove_non_intersect_data(data,dir_i,ind_intersect);


%solo bsp
ind_data={ind_intersect.solo_bsp_ind_for_analysis};
x=[ind_intersect.solo_bsp_ind_for_analysis];

non_empty=find(~cellfun(@isempty,ind_data));
size_mat=size(ind_data{non_empty(1)});
all_mats=reshape(x,size_mat(1),size_mat(2),[]);
logical_ind=sum(all_mats,3);
ind_to_remove=~logical_ind;

data(dir_i).solo.bsp.flight_ts(ind_to_remove)=nan;
data(dir_i).solo.bsp.flight_x_pos(ind_to_remove)=nan;
data(dir_i).solo.bsp.flight_y_pos(ind_to_remove)=nan;

%solo spikes
ind_data={ind_intersect.solo_spike_ind_for_analysis};
x=[ind_intersect.solo_spike_ind_for_analysis];

non_empty=find(~cellfun(@isempty,ind_data));
size_mat=size(ind_data{non_empty(1)});
all_mats=reshape(x,size_mat(1),size_mat(2),[]);
logical_ind=sum(all_mats,3);
ind_to_remove=~logical_ind;

data(dir_i).solo.spikes.flight_ts(ind_to_remove)=nan;
data(dir_i).solo.spikes.flight_x_pos(ind_to_remove)=nan;
data(dir_i).solo.spikes.flight_y_pos(ind_to_remove)=nan;

%co bsp
ind_data={ind_intersect.co_bsp_ind_for_analysis};
x=[ind_intersect.co_bsp_ind_for_analysis];

non_empty=find(~cellfun(@isempty,ind_data));
size_mat=size(ind_data{non_empty(1)});
all_mats=reshape(x,size_mat(1),size_mat(2),[]);
logical_ind=sum(all_mats,3);
ind_to_remove=~logical_ind;

data(dir_i).co.bsp.flight_ts(ind_to_remove)=nan;
data(dir_i).co.bsp.flight_x_pos(ind_to_remove)=nan;
data(dir_i).co.bsp.flight_y_pos(ind_to_remove)=nan;

%co spikes
ind_data={ind_intersect.co_spike_ind_for_analysis};
x=[ind_intersect.co_spike_ind_for_analysis];

non_empty=find(~cellfun(@isempty,ind_data));
if ~isempty(non_empty)
size_mat=size(ind_data{non_empty(1)});
all_mats=reshape(x,size_mat(1),size_mat(2),[]);
logical_ind=sum(all_mats,3);
ind_to_remove=~logical_ind;

data(dir_i).co.spikes.flight_ts(ind_to_remove)=nan;
data(dir_i).co.spikes.flight_x_pos(ind_to_remove)=nan;
data(dir_i).co.spikes.flight_y_pos(ind_to_remove)=nan;
data(dir_i).co.spikes.flight_dis(ind_to_remove)=nan;
end
