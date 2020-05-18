function    [shuffle_solo_bsp,shuffle_solo_spikes]=PV_find_solo_data_for_shuffle(solo_data,co_data,ego_bin_start,ego_bin_end,bin_i,dir_i)

%get relevant co data:
co_dis_bsp=co_data.bsp.flight_dis';co_dis_bsp=co_dis_bsp(:);
co_x_pos=co_data.bsp.flight_x_pos';co_x_pos=co_x_pos(:)';

relevant_ind_bsp=(co_dis_bsp>=ego_bin_start(bin_i) & co_dis_bsp<ego_bin_end(bin_i));

%find co epochs
[ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(find(relevant_ind_bsp),length(relevant_ind_bsp));

bsp_data_flight_to_remove_ind_all=[];
spike_data_flight_to_remove_ind_all=[];

for epoch_i=1:length(ind_length)% run over all CO that had data in this distance bin
    
    %find the epoch positions
    epoch_pos=co_x_pos(start_ind(epoch_i):end_ind(epoch_i));
    % find solo data flight with this positions
    if dir_i==1
        [row,col]=find(solo_data.bsp.flight_x_pos>=epoch_pos(1) & solo_data.bsp.flight_x_pos<=epoch_pos(end));
    else
        [row,col]=find(solo_data.bsp.flight_x_pos>=epoch_pos(end) & solo_data.bsp.flight_x_pos<=epoch_pos(1));
        
    end
    if isempty(row)
        continue
    end
    relevanlt_flights=unique(row);
    %choose random flight:
    flight_to_remove_data_ind=randi(length(relevanlt_flights));
    flight_chosen=relevanlt_flights(flight_to_remove_data_ind);
    %take the data to shuffle i bsp:
    len_epoch_to_remove=length(col(row==flight_chosen));
    bsp_data_flight_to_remove_ind=[repmat(flight_chosen,len_epoch_to_remove,1),col(row==flight_chosen)];
    %take the data to shuffle i spikes:
    spikes_pos_flight_i=solo_data.spikes.flight_x_pos(flight_chosen,:);
    if dir_i==1
        
        ind_spike_to_remove=find(spikes_pos_flight_i>=epoch_pos(1) & spikes_pos_flight_i<=epoch_pos(end));
    else
        ind_spike_to_remove=find(spikes_pos_flight_i>=epoch_pos(end) & spikes_pos_flight_i<=epoch_pos(1));
        
    end
    if ~isempty(ind_spike_to_remove)
        spike_data_flight_to_remove_ind=[repmat(flight_chosen,length(ind_spike_to_remove),1),ind_spike_to_remove'];
    else
        spike_data_flight_to_remove_ind=[];
    end
    bsp_data_flight_to_remove_ind_all=[bsp_data_flight_to_remove_ind_all;bsp_data_flight_to_remove_ind];
    spike_data_flight_to_remove_ind_all=[spike_data_flight_to_remove_ind_all;spike_data_flight_to_remove_ind];
    
end

if ~isempty(bsp_data_flight_to_remove_ind_all)
shuffle_solo_bsp=solo_data.bsp.flight_x_pos(sub2ind(size(solo_data.bsp.flight_x_pos),bsp_data_flight_to_remove_ind_all(:,1),bsp_data_flight_to_remove_ind_all(:,2)));
else
   shuffle_solo_bsp=[]; 
end
if ~isempty(spike_data_flight_to_remove_ind_all)
    shuffle_solo_spikes=solo_data.spikes.flight_x_pos(sub2ind(size(solo_data.spikes.flight_x_pos),spike_data_flight_to_remove_ind_all(:,1),spike_data_flight_to_remove_ind_all(:,2)));
else
    shuffle_solo_spikes=[];
    
end