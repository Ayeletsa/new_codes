function calc_inclusion_criteria(inclusion_cells_params_param_file_name)

%% load params
load(inclusion_cells_params_param_file_name)


files = dir(cell_co_solo_initial_analysis_struct_folder);
cell_co_solo_initial_struct_names = {files.name};
cell_co_solo_initial_struct_names([files.isdir])=[];
ego_signif=[];
place_cells_signif=[];
valid_cells=[];
%% for each cell
for ii_cell = 1:length(cell_co_solo_initial_struct_names)
    ii_cell
    % load cell's egocentric struct
    struct_name = cell_co_solo_initial_struct_names{ii_cell};
    file_name = fullfile(cell_co_solo_initial_analysis_struct_folder,struct_name);
    load(file_name);
    
    % save basic information about the cell
    bat = cell_co_solo_initial_analysis.exp_data.bat;
    day = cell_co_solo_initial_analysis.exp_data.day;
    if isnumeric(day)
        day=num2str(day);
    end
    cell_num = cell_co_solo_initial_analysis.exp_data.cell_num;
    
    % shuffle names:
    solo_shuffle_file_name= fullfile(solo_shuffle_folder_name,['solo_shuffling_struct_b',num2str(bat),'_d', day,'_c',num2str(cell_num)]);
    co_shuffle_file_name= fullfile(co_shuffle_folder_name,['co_shuffling_struct_b',num2str(bat),'_d', day,'_c',num2str(cell_num)]);
    inclusion=struct();
    for ii_dir = 1:2
        %% 1. check if a cell is valid
        exp_data=cell_co_solo_initial_analysis.exp_data;
        inclusion=check_if_valid(inclusion,exp_data,ii_dir,inclusion_cells_params_param_file_name);
        %% 2. check if signif place cell:
        solo=cell_co_solo_initial_analysis.solo;
        load(solo_shuffle_file_name);
        inclusion=check_if_place_cell(inclusion,solo,ii_dir,shuffling_struct,inclusion_cells_params_param_file_name);
        %% 3. check if signif egocentric:
        co=cell_co_solo_initial_analysis.co;
        load(co_shuffle_file_name);
        inclusion=check_if_ego_cell(inclusion,co,ii_dir,shuffling_struct,inclusion_cells_params_param_file_name);
    end
    
    
    
    %% save data to file
    filename = fullfile(inclusion_cells_folder_name,['inclusion_cell_',num2str(cell_num)]);
    save(filename, 'inclusion');
    
    %% add to all cells inclusion:
    if sum([inclusion.ego_cell])>0
        ego_signif=[ego_signif,cell_num];
    end
    if sum([inclusion.place_cell])>0
        place_cells_signif=[place_cells_signif,cell_num];
    end
    if sum([inclusion.valid_cell])>0
        valid_cells=[valid_cells,cell_num];
    end
    
end
%% save all inclusion:
all_cell_inclusion_list.ego_signif=ego_signif;
all_cell_inclusion_list.place_cells_signif=place_cells_signif;
all_cell_inclusion_list.valid_cells=valid_cells;

filename = fullfile(inclusion_cells_folder_name,['all_cell_inclusion_list']);
save(filename, 'all_cell_inclusion_list');

end