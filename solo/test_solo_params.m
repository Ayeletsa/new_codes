clear
% define parameters:
kernels=[0.5,1,1.5];
for ii=1:3
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
kernel=kernels(ii);
Twobatproj_param(param_folder,kernel)


behav_param_file_name=fullfile(param_folder,'behav_params.mat');
dir_param_file_name=fullfile(param_folder,'dirs_params.mat');
solo_param_file_name=fullfile(param_folder,'solo_params.mat');
co_param_file_name=fullfile(param_folder,'co_params.mat');
population_param_file_name=fullfile(param_folder,'co_population_params.mat');
co_shuffle_param_file_name=fullfile(param_folder,'co_shuffle_params.mat');

field_param_file_name=fullfile(param_folder,'fields_params.mat');
%% 2. initial CO and solo anlysis
initial_co_solo_analysis(behav_param_file_name,dir_param_file_name,solo_param_file_name,co_param_file_name,field_param_file_name)


end