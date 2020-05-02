clear
% define parameters:
kernels=[0.5,1,1.5];
FR_vec=[0.5,1];
for bat=1:2
    if bat==1
        cells_struct_dir='L:\Data\2batproj\Data_Nlg_Proc\yr_2019_bat_2336\cell_structs\';
        
        
    else
        cells_struct_dir='D:\Ayelet\Data\Data_Nlg_Proc\yr_2018_bat_2389\cell_structs\';
    end
    for ki=1:3
        for fr_i=1:2
            FR_thr=(FR_vec(fr_i));
            param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
            kernel=kernels(ki);
            Twobatproj_param_test_solo_params(param_folder,kernel,cells_struct_dir,FR_thr)
            
            
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
    end
end