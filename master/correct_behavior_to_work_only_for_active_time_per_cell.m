function [logical_vec_of_active_flight]=correct_behavior_to_work_only_for_active_time_per_cell(cell_struct,bsp_ts_usec)


%% corect behavior fo relevant time for cell
relevant_timestamps_for_cell=cell_struct.cell_info.relevant_timestamps_for_cell;
behave_session=2;
behave_ts=[cell_struct.exp_info.nlg_events(behave_session).start_time cell_struct.exp_info.nlg_events(behave_session).end_time];
ind=1;
corect_behave=0;
if relevant_timestamps_for_cell(ind)>behave_ts(ind)
    cell_ts(ind)=relevant_timestamps_for_cell(ind);
    corect_behave=1;
    correct_from_begin=1;
    
else
    cell_ts(ind)=behave_ts(ind);
end
ind=2;
if relevant_timestamps_for_cell(ind)<behave_ts(ind)
    cell_ts(ind)=relevant_timestamps_for_cell(ind);
    corect_behave=1;
    correct_from_end=1;
    
else
    cell_ts(ind)=behave_ts(ind);
end
%%
logical_vec_of_active_flight=ones(length(bsp_ts_usec),1);
if corect_behave==1
  % find flight where cell was not active:
 
    flight_ind_to_remove= find(~cellfun(@isempty,cellfun(@(x) (find(x<=cell_ts(1) | x>=cell_ts(end))), bsp_ts_usec, 'UniformOutput', false)));
   logical_vec_of_active_flight(flight_ind_to_remove)=0;
end