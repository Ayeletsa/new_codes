function initial_co_solo_analysis(behav_param_file_name,dir_param_file_name,solo_param_file_name,co_param_file_name,field_param_file_name,per_field_param_file_name)
load(behav_param_file_name)
load(dir_param_file_name)
data_dir_info=dir(cells_struct_dir);


%% for each cell


for ii_cell = 3:length(data_dir_info)
    
    
    ii_cell
    % load cell struct
    %% load cell struct:
    file_name=fullfile(cells_struct_dir,data_dir_info(ii_cell).name)
    load(file_name)
    
    % get basic information
    bat = cell_struct.exp_info.bat;
    day = cell_struct.cell_info.day;
    tag_self = cell_struct.exp_info.bsp_tag_self;
    tags = [cell_struct.bsp_data.tag_ID];
    tag_i = find(ismember(tags,tag_self));
    behav_struct_name=[behave_day_struct_folder,'behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat'];
    general_behavior_data_file_name=fullfile(behave_day_struct_folder,['general_behave_analysis_bat',num2str(bat),'_day',num2str(day),'.mat']);

    % load behvioral mode struct
    load(behav_struct_name)
    load(general_behavior_data_file_name)
    
 
    %%
    % save into behavior struct
    cell_co_solo_initial_analysis.exp_data.bat = bat;
    cell_co_solo_initial_analysis.exp_data.day = day;
    cell_co_solo_initial_analysis.exp_data.cell_num = cell_struct.cell_info.cell_num  ;
    [cell_co_solo_initial_analysis.exp_data.stability, cell_co_solo_initial_analysis.exp_data.mean_fr,cell_co_solo_initial_analysis.exp_data.n_spike_per_session] = stability_index(cell_struct);
    cell_co_solo_initial_analysis.exp_data.spikes_ts_usec=cell_struct.spikes.spikes_ts_usec;
    cell_co_solo_initial_analysis.exp_data.L_Ratio=cell_struct.spikes.L_Ratio;
    cell_co_solo_initial_analysis.exp_data.Isolation_dis=cell_struct.spikes.Isolation_dis;
    % analyze each behavior separately
    solo_struct_name=fullfile(behave_solo_struct_folder,['solo_bsp_data_bat_',num2str(bat),'_day_',num2str(day)]);
    cell_co_solo_initial_analysis.solo = spike_struct_solo_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_param_file_name,field_param_file_name,solo_struct_name);
    
    co_struct_name=fullfile(behave_co_struct_folder,['co_bsp_data_bat_',num2str(bat),'_day_',num2str(day)]);
    cell_co_solo_initial_analysis.co = spike_struct_co_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i, cell_co_solo_initial_analysis.solo,co_param_file_name,per_field_param_file_name,field_param_file_name,bat,day,co_struct_name);
    %behavior_struct.obs = behavior_struct_obs_data(bsp_data,cell_struct,tag_i);
    %behavior_struct.tr = behavior_struct_tracking_data(bsp_data,cell_struct,behavioral_modes,tag_i);
    
    file_name = [cell_co_solo_initial_analysis_struct_folder,'\bat',num2str(bat),'_day_',num2str(day),'_cell_', num2str(cell_struct.cell_info.cell_num),'.mat'];
    
    save(file_name,'cell_co_solo_initial_analysis')
    clear behavior_struct cell_co_solo_initial_analysis
    
end