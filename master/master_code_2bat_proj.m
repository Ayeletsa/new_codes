%% Master code 2d proj
clear
% define parameters:
[cells_struct_dir,main_analysis_dir,param_folder]=local_dirs;
Twobatproj_param(cells_struct_dir,main_analysis_dir,param_folder)

behav_param_file_name=fullfile(param_folder,'behav_params.mat');
dir_param_file_name=fullfile(param_folder,'dirs_params.mat');
solo_param_file_name=fullfile(param_folder,'solo_params.mat');
co_param_file_name=fullfile(param_folder,'co_params.mat');
population_param_file_name=fullfile(param_folder,'co_population_params.mat');
co_shuffle_param_file_name=fullfile(param_folder,'co_shuffle_params.mat');
solo_shuffle_param_file_name=fullfile(param_folder,'solo_shuffle_params.mat');
field_param_file_name=fullfile(param_folder,'fields_params.mat');
per_field_param_file_name=fullfile(param_folder,'per_field_params.mat');
population_vector_param_file_name=fullfile(param_folder,'population_vector_params.mat');
inclusion_cells_params_param_file_name=fullfile(param_folder,'inclusion_cells_params.mat');

%% 1. analyze behavior 
create_behavior_structs(behav_param_file_name,dir_param_file_name)
%% 2. initial CO and solo anlysis
initial_co_solo_analysis(behav_param_file_name,dir_param_file_name,solo_param_file_name,co_param_file_name,field_param_file_name,per_field_param_file_name)
%% 3. CO shuffling
co_shuffling_new(co_shuffle_param_file_name)
%% 4. Solo shuffling
solo_shuffling(solo_shuffle_param_file_name,field_param_file_name)
%% 5. CO figures
plot_co_main_cell_fig(co_param_file_name,dir_param_file_name,population_param_file_name,per_field_param_file_name)
%% 6. co figure per field time
plot_co_main_cell_fig_time_per_field(dir_param_file_name,population_param_file_name,co_param_file_name,per_field_param_file_name)

%% 7.CO figure with per inter field analysis
plot_co_main_cell_fig_new_inter_field(dir_param_file_name,population_param_file_name,co_param_file_name,per_field_param_file_name)


%% create inclustion lists:
calc_inclusion_criteria(inclusion_cells_params_param_file_name)