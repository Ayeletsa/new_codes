function create_behavior_structs(behav_param_file_name,dir_param_file_name,solo_param_file_name,co_param_file_name)

load(behav_param_file_name)
load(dir_param_file_name)
data_dir_info=dir(cells_struct_dir);


%% for each cell
day_before=0;

for ii_cell = 542:length(data_dir_info)

    
    ii_cell
    
    %% load cell struct:
    file_name=fullfile(cells_struct_dir,data_dir_info(ii_cell).name);
    load(file_name)
    
    % get basic information
    bat = cell_struct.exp_info.bat;
    %day = cell_struct.exp_info.day;
    day =cell_struct.cell_info.day;
    if isnumeric(day)
        day=num2str(day);
        
    end
    tag_self = cell_struct.exp_info.bsp_tag_self;
    tags = [cell_struct.bsp_data.tag_ID];
    tag_i = find(ismember(tags,tag_self));
    
    behav_struct_name=[behave_day_struct_folder,'behavioral_modes_bat_',num2str(bat),'_day_',day,'.mat'];
    
    % analyze behavior
    behave_ts=[cell_struct.exp_info.nlg_events(2).start_time,cell_struct.exp_info.nlg_events(2).end_time];
    
    ball_pos_name=[ball_position_folder,'ball_pos_bat_',num2str(bat),'_day_',day,'.mat'];
    bsp_data=cell_struct.bsp_data;
       %start_end_ts_ns = [cell_struct.exp_info.nlg_events(2).start_time_fitted_to_bsp_msec cell_struct.exp_info.nlg_events(3).end_time_fitted_to_bsp_msec] * 1e6;
    try
        load(behav_struct_name) 
    catch
        if ~strcmp(day,day_before)
            [behavioral_modes]=find_behavioral_modes(bsp_data ,behav_param_file_name,tag_i,bat,day,behave_ts,dir_param_file_name,ball_pos_name); % assign bsp samples to different behaviors
        else
            load(behav_struct_name)
        end
   
    end
    
    clear behavior_struct
    day_before=day;
end

end