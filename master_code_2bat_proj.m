%% Master code 2d proj
clear
% define parameters:
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
Twobatproj_param(param_folder)

%% 1. analyze behavior 
behav_param_file_name=fullfile(param_folder,'behav_params.mat');
dir_param_file_name=fullfile(param_folder,'dirs_params.mat');
solo_param_file_name=fullfile(param_folder,'solo_params.mat');
co_param_file_name=fullfile(param_folder,'co_params.mat');

create_behavior_structs(behav_param_file_name,dir_param_file_name)
%% 2. initial CO and solo anlysis
initial_co_solo_analysis(behav_param_file_name,dir_param_file_name,solo_param_file_name,co_param_file_name)
%% 3. CO shuffling
co_shuffle_param_file_name=fullfile(param_folder,'co_shuffle_params.mat');

co_shuffling_new(co_shuffle_param_file_name)

%% 4. Solo shuffling


%% 5. CO figures
population_param_file_name=fullfile(param_folder,'co_population_params.mat');

plot_co_main_cell_fig(dir_param_file_name,population_param_file_name)
