function clicks_spikes_analysis(dir_param_file_name,audio_param_file_name,co_param_file_name,solo_param_file_name)
load(co_param_file_name)
load(dir_param_file_name)
data_dir_info=dir(cells_struct_dir);
spikes_files = dir(cell_co_solo_initial_analysis_struct_folder);
shuffle_files = dir(co_shuffle_folder_name);

%% for each day
previous_day = '';
for ii_cell = 3:length(data_dir_info)
    ii_cell
    
    %% load cell struct:
    file_name=fullfile(cells_struct_dir,data_dir_info(ii_cell).name);
    load(file_name)
    [~,cell_name,~] = fileparts(file_name);
    
    % use only cells which have audio
    if ~cell_struct.exp_info.audio_self
        continue
    end
    
    % get basic information
    bat = cell_struct.exp_info.bat;
    day = cell_struct.exp_info.day;
    name_prefix = ['bat_',num2str(bat),'_day_',num2str(day)];
    tag_self = cell_struct.exp_info.bsp_tag_self;
    tags = [cell_struct.bsp_data.tag_ID];
    tag_i = find(ismember(tags,tag_self));
    behav_struct_name=[behave_day_struct_folder,'behavioral_modes_',name_prefix,'.mat'];
    general_behavior_data_file_name=fullfile(behave_day_struct_folder,['general_behave_analysis_bat',num2str(bat),'_day',num2str(day),'.mat']);
    
    % load behvioral mode struct
    load(behav_struct_name)
    load(general_behavior_data_file_name,'bsp_proc_data')
    % save into behavior struct
    cell_clicks_spikes_analysis.exp_data.bat = bat;
    cell_clicks_spikes_analysis.exp_data.day = day;
    cell_clicks_spikes_analysis.exp_data.cell_num = cell_struct.cell_info.cell_num;
    [cell_clicks_spikes_analysis.exp_data.stability, cell_clicks_spikes_analysis.exp_data.mean_fr] = stability_index(cell_struct);
    
    %% create clicks struct (once for each day, so no need to re-do it for every cell) 
    if ~strcmp(num2str(previous_day),num2str(day))
        clicks_struct_name=[clicks_day_struct_folder,'clicks_',name_prefix,'.mat'];
        audio_filt_name = [audio_struct_dir,name_prefix,'_audio_filt'];
        create_clicks_day_struct(clicks_struct_name,audio_filt_name,audio_param_file_name,bsp_proc_data,tag_i,behavioral_modes,co_param_file_name,solo_param_file_name,click_detection_folder);
        %         clicks_struct(1).co = clicks_in_co(bsp_proc_data,behavioral_modes,tag_i,co_param_file_name,clicks_struct,0);
        %         clicks_struct(1).solo = clicks_in_solo(bsp_proc_data,behavioral_modes,tag_i,solo_param_file_name,clicks_struct,0);
        %         save(clicks_struct_name,'clicks_struct');
        
        %% divide COs by the calling rate increase latency
%         cr = divide_early_late(clicks_struct_name);
%         cell_clicks_spikes_analysis.calling_rate = cr;
    else
%         cell_clicks_spikes_analysis.calling_rate = cr;
    end
    
    previous_day = day;
    continue
    
    %% analyze spikes+clicks
    % load cell structs
    load(fullfile(cell_co_solo_initial_analysis_struct_folder,spikes_files(ii_cell).name));
    load(fullfile(co_shuffle_folder_name,shuffle_files(ii_cell).name));
    bsp = cell_co_solo_initial_analysis.co(ii_dir).bsp;
    spikes = cell_co_solo_initial_analysis.co(ii_dir).spikes;
    % decide on criterion!!!
    prob = problem_with_solo
    thresh_used = thresh_solo;
    co_early = rate_inc_solo < thresh_used;
    

    
    %     bsp = structfun(@(x) x(~high_bl & ~prob,:),bsp,'UniformOutput',false);
    %     spikes = structfun(@(x) x(~high_bl & ~prob,:),spikes,'UniformOutput',false);

            
            
            figure
            clf
            p=panel();
            c=corr(dis_x_pos_fr_for_2D_early',dis_x_pos_fr_for_2D_late');
            sgtitle([cell_name ', dir_' num2str(ii_dir) ', corr=' num2str(c)],'interpreter','none')
            p.pack({0.2,[]},2);
            %     p.select('all')
            %     p.identify();
            p(1,1).select()
            title('early')
            plot(firing_rate.dis_X_bins_vector_of_centers{1,1},dis_x_pos_fr_for_2D_early,'m','LineWidth',2)
            p(2,1).select()
            fn_plot_2D_field(allo_ego_map_early,firing_rate(ii_dir).dis_X_bins_vector{1},firing_rate(ii_dir).dis_X_bins_vector_of_centers{1},firing_rate(ii_dir).allo_X_bins_vector{1},firing_rate(ii_dir).allo_X_bins_vector_of_centers{1},[]);
            axis xy
            p(1,2).select()
            title('late')
            plot(firing_rate.dis_X_bins_vector_of_centers{1,1},dis_x_pos_fr_for_2D_late,'m','LineWidth',2)
            p(2,2).select()
            fn_plot_2D_field(allo_ego_map_late,firing_rate(ii_dir).dis_X_bins_vector{1},firing_rate(ii_dir).dis_X_bins_vector_of_centers{1},firing_rate(ii_dir).allo_X_bins_vector{1},firing_rate(ii_dir).allo_X_bins_vector_of_centers{1},[]);
            axis xy
            p.margin = [15 15 15 20];
            
            %     bsp_mat_allo = bsp.x_pos(high_bl,:);
            %     bsp_vec_allo = bsp_mat_allo(isfinite(bsp_mat_allo));
            %     spikes_mat_allo = spikes.x_pos(high_bl,:);
            %     spikes_vec_allo = spikes_mat_allo(isfinite(spikes_mat_allo));
            %     bsp_mat_dis = bsp.dis_m(high_bl,:);
            %     bsp_vec_dis = bsp_mat_dis(isfinite(bsp_mat_dis));
            %     spikes_mat_dis = spikes.dis_m(high_bl,:);
            %     spikes_vec_dis = spikes_mat_dis(isfinite(spikes_mat_dis));
            %
            %     [~, ~, allo_ego_map_high_bl, ~, ~, ~, ~, ~, ~, ~] ...
            %             = fn_compute_generic_2D_field ...
            %             (firing_rate.dis_X_bins_vector{1,1}, firing_rate.allo_X_bins_vector{1,1}, bsp_vec_dis, bsp_vec_allo, spikes_vec_dis, spikes_vec_allo, time_spent_minimum_for_2D_bins, frames_per_second, sigma_a, hsize, legalize_by_neighbor_bins_flag);
            %
            % analyze each behavior separately
            %     cell_clicks_spikes_analysis.solo = spike_struct_solo_data (bsp_proc_data,cell_struct,behavioral_modes,tag_i,solo_param_file_name);
            %
            %     file_name = [clicks_spikes_analysis_struct_folder,'\bat',num2str(bat),'_day_',num2str(day),'_cell_', num2str(cell_struct.cell_info.cell_num),'.mat'];
            %
            %     save(file_name,'clicks_spikes_analysis')
            
end
