function shuffles = rigid_shuffle_of_trials(bsp,spikes,nshuffles,min_offset_perc)
sp_mat_transpose = spikes.ts_usec';
b_mat_transpose = bsp.ts_usec';
spikes_ts = sp_mat_transpose(isfinite(sp_mat_transpose));
behav_ts = b_mat_transpose(isfinite(b_mat_transpose));
behav_ind = find(isfinite(b_mat_transpose));

fields = fieldnames(spikes);
min_offset=round(min_offset_perc*length(behav_ts));
p = randi([min_offset,length(behav_ts)-min_offset],nshuffles,1);
p=[0;p];
circular_behav_ind = repmat(behav_ind,2,1);
shuffles = struct();

%spikes in behavior ind
[~,spikes_ind] = arrayfun(@(x) min(abs(x-behav_ts)),spikes_ts);

parfor ishuffle = 1:nshuffles+1
%     if ishuffle==1
%         shift(ishuffle)=0;
%     else
%         shift(ishuffle)=p(ishuffle-1);
%     end
    new_spikes_ind_circ = circular_behav_ind(spikes_ind+p(ishuffle));
    for ifield=1:length(fields)
        cur_field = bsp.(fields{ifield})';
        shuffles(ishuffle).(fields{ifield}) = cur_field(new_spikes_ind_circ);
    end
end

