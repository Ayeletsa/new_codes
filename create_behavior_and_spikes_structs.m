function create_behavior_and_spikes_structs(behav_param_file_name,dir_param_file_name,solo_param_file_name,co_param_file_name)

load(behav_param_file_name)
load(dir_param_file_name)
data_dir_info=dir(cells_struct_dir);


%% for each cell
day_before=0;

for ii_cell = 70:length(data_dir_info)

    
    ii_cell
    
    %% load cell struct:
    file_name=fullfile(cells_struct_dir,data_dir_info(ii_cell).name);
    load(file_name)
    
    % get basic information
    bat = cell_struct.exp_info.bat;
    day = cell_struct.exp_info.day;
    tag_self = cell_struct.exp_info.bsp_tag_self;
    tags = [cell_struct.bsp_data.tag_ID];
    tag_i = find(ismember(tags,tag_self));
    
    behav_struct_name=[behave_day_struct_folder,'behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat'];
    
    % analyze behavior
    behave_ts=[cell_struct.exp_info.nlg_events(2).start_time,cell_struct.exp_info.nlg_events(2).end_time];
    ball_pos_name=[ball_position_folder,'ball_pos_bat_',num2str(bat),'_day_',num2str(day),'.mat'];
    
    % interpolate position of the other bat to have the same ts for both
    % bats
     ts=cell_struct.bsp_data(tag_i).ts_us_upsampled;
     ts_original=cell_struct.bsp_data(3-tag_i).ts_us_upsampled;
     pos_other_original=cell_struct.bsp_data(3-tag_i).pos_upsampled;
     pos_other_x=interp1(ts_original,pos_other_original(1,:),ts);
     pos_other_y=interp1(ts_original,pos_other_original(2,:),ts);

     %define relevant bsp_data struct:
    bsp_data(3-tag_i).pos=[pos_other_x;pos_other_y]';
    bsp_data(3-tag_i).ts=ts';
    
    bsp_data(tag_i).pos=cell_struct.bsp_data(tag_i).pos_upsampled';
    bsp_data(tag_i).ts=cell_struct.bsp_data(tag_i).ts_us_upsampled';

    bsp_data=find_flight_ind(bsp_data,behave_ts,behav_param_file_name,ball_pos_name); % find when the bat is flying
    %start_end_ts_ns = [cell_struct.exp_info.nlg_events(2).start_time_fitted_to_bsp_msec cell_struct.exp_info.nlg_events(3).end_time_fitted_to_bsp_msec] * 1e6;
    try
        load(behav_struct_name) 
    catch
        if day~=day_before
            [behavioral_modes]=find_behavioral_modes(bsp_data ,behav_param_file_name,tag_i,bat,day,behave_ts,dir_param_file_name); % assign bsp samples to different behaviors
        else
            load(behav_struct_name)
        end
   
    end
    %load(behav_struct_name)
    % save into behavior struct
    behavior_struct.exp_data.bat = bat;
    behavior_struct.exp_data.day = day;
    behavior_struct.exp_data.cell_num = cell_struct.cell_info.cell_num  ;
    [behavior_struct.exp_data.stability, behavior_struct.exp_data.mean_fr] = stability_index(cell_struct);
    
    % analyze each behavior separately
    behavior_struct.solo = spike_struct_solo_data (bsp_data,cell_struct,behavioral_modes,tag_i,solo_param_file_name);
    behavior_struct.co = spike_struct_co_data (bsp_data,cell_struct,behavioral_modes,tag_i, behavior_struct.solo,co_param_file_name);
    %behavior_struct.obs = behavior_struct_obs_data(bsp_data,cell_struct,tag_i);
    %behavior_struct.tr = behavior_struct_tracking_data(bsp_data,cell_struct,behavioral_modes,tag_i);
    
    file_name = [behave_cell_struct_folder,'\bat',num2str(bat),'_day_',num2str(day),'_cell_', num2str(cell_struct.cell_info.cell_num),'.mat'];
    
    save(file_name,'behavior_struct')
    clear behavior_struct
    day_before=day;
end

end