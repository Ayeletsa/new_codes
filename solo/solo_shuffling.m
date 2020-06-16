function solo_shuffling(solo_shuffle_param_file_name,field_param_file_name)

load(solo_shuffle_param_file_name)
prm.fields=load(field_param_file_name);


files = dir(cell_co_solo_initial_analysis_struct_folder);
behavior_struct_names = {files.name};
%% for each cell
for ii_cell = 3:length(behavior_struct_names)
    ii_cell
    % load cell's egocentric struct
    struct_name = behavior_struct_names{ii_cell};
    file_name = fullfile(cell_co_solo_initial_analysis_struct_folder,struct_name);
    load(file_name);
    
    % save basic information about the cell
    shuffling_struct(1).info.bat = cell_co_solo_initial_analysis.exp_data.bat;
    shuffling_struct(1).info.day = cell_co_solo_initial_analysis.exp_data.day;
    if isnumeric(shuffling_struct(1).info.day)
        shuffling_struct(1).info.day=num2str(shuffling_struct(1).info.day);      
    end
    shuffling_struct(1).info.cell_num = cell_co_solo_initial_analysis.exp_data.cell_num;
    
    % shuffle each fircetion separately
    for ii_dir = 1:2
        FE=cell_co_solo_initial_analysis.solo(ii_dir).FE_struct_for_tamirs_code;
        if ~isempty(FE)
        shuffling_struct(ii_dir).FE_PSTH_shuffle = FE_compute_PSTH_shuffle(FE,...
                                                      shuffles_num,...
                                                      shuffles_max_shift,prm);
        else
            shuffling_struct(ii_dir).FE_PSTH_shuffle =nan;
        end
    end
    
    %% save data to structs:
    shuffle_struct_name = ['solo_shuffling_struct_b',num2str( shuffling_struct(1).info.bat),'_d', shuffling_struct(1).info.day,'_c',num2str(shuffling_struct(1).info.cell_num)];
    file_name=fullfile(solo_shuffle_folder_name,shuffle_struct_name);
    save(file_name,'shuffling_struct')
    clear shuffling_struct
end