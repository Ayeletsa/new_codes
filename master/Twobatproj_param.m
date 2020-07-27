function Twobatproj_param(cells_struct_dir,main_analysis_dir,param_folder)

%% folders:
%load('local_computer_dirs')
% % data input:
 params.dirs.cells_struct_dir=cells_struct_dir;
% %params.dirs.cells_struct_dir='L:\Data\2batproj\Data_Nlg_Proc\yr_2019_bat_2336\cell_structs\';
% %params.dirs.cells_struct_dir='D:\Ayelet\Data\Data_Nlg_Proc\yr_2018_bat_2389\cell_structs\';
% %params.dirs.cells_struct_dir='L:\Data\2batproj\Data_Nlg_Proc\yr_2019_bat_2299';
% 
% % data output:
% main_analysis_dir='D:\Ayelet\2bat_proj\Analysis\new_code\';
params.dirs.behave_day_struct_folder=[main_analysis_dir,'\analysis_structs\behavioral_modes\day_structs\'];
params.dirs.behave_co_struct_folder=[main_analysis_dir,'\analysis_structs\behavioral_modes\day_structs\co_structs\'];
params.dirs.behave_solo_struct_folder=[main_analysis_dir,'\analysis_structs\behavioral_modes\day_structs\solo_structs\'];
params.dirs.ball_position_folder=[main_analysis_dir,'\analysis_structs\behavioral_modes\ball_pos\'];
params.dirs.cell_co_solo_initial_analysis_struct_folder=[main_analysis_dir,'\analysis_structs\co_solo_initial_analysis\'];
params.dirs.solo_initial_analysis_struct_folder=[main_analysis_dir,'\analysis_structs\solo_initial\'];
params.dirs.co_initial_analysis_struct_folder=[main_analysis_dir,'\analysis_structs\co_initial\'];
%params.dirs.cell_co_solo_initial_analysis_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis_k_1.5_th_1\';
params.dirs.co_shuffle_folder_name = [main_analysis_dir,'\analysis_structs\co_shuffling_struct'];
params.dirs.solo_shuffle_folder_name=[main_analysis_dir,'\analysis_structs\solo_shuffling_struct'];
params.dirs.inclusion_cells_folder_name=[main_analysis_dir,'\analysis_structs\inclusion_cells_struct'];
% figures:
params.dirs.behave_analysis_fig_dir_out=[main_analysis_dir,'figures\initial_behavior_analysis\'];
params.dirs.co_fig_folder_name=[main_analysis_dir,'figures\basic_co_analysis'];
params.dirs.co_signif_cells_fig_folder_name=[params.dirs.co_fig_folder_name,'signif_cells\'];

% create folders if does not exist
dirs=fieldnames(params.dirs);
for dir_i=1:length(dirs)
  if ~exist(params.dirs.(dirs{dir_i}))
      mkdir(params.dirs.(dirs{dir_i}));
  end
end

if ~exist(param_folder)
    mkdir(param_folder)
end
%save
dir_params=params.dirs;
param_file_name=fullfile(param_folder,'dirs_params.mat');
save(param_file_name, '-struct', 'dir_params')

%% parameters for find_flight_ind
params.behav.min_velocity_flight=2; % define flight by velocity (m/sec)
params.behav.dist_from_the_ball=3;
params.behav.min_time_for_slow_flights=1*100;
% ball position - later you need to find the exact ball position for each
% day, but this is the area to look for the ball. make sure to change if
% needed.
params.behav.ball_1_area=[0 10];
params.behav.ball_2_area=[110 140];
params.behav.landing_vel=0.2;%m/s
params.behav.flight_upper_vel_thres=30; %m/s
%% parametres for cross-overs analysis
params.behav.time_before_after_co = 3; % in secs
params.behav.dis_before_after_co = 40; % in meters

%% parameters for behavioral modes
params.behav.load_struct=1; %1 if you want to load behavioral struct
params.behav.correct_manually=1;
params.behav.frame_per_second=100;
params.behav.csaps_p = 1e-5; %for vel computation

params.behav.min_solo_length=2*params.behav.frame_per_second; %samples
params.behav.min_tracking_length=3*params.behav.frame_per_second; %samples
params.behav.dist_thresh_tracking=20; %meters
params.behav.dist_thresh_solo=40; %meters
params.behav.min_dist_opposite_dirs_before_after_CO=5; %meters
params.behav.time_before_after_co_for_co_window=3.5; %(seconds) for defining the window around the co
% for computing reduction in velocity:
params.behav.long_dis_thresh=20;
params.behav.short_dis_thresh=5;

%params.behav.CO_window=[20 20]; %meters
%params.behav.max_wind_CO=2*params.behav.frame_per_second;
%params.behav.min_time_before_CO=2*params.behav.frame_per_second;
%params.behav.UT_window=params.behav.frame_per_second;
%params.behav.UT_time_from_CO=2*params.behav.frame_per_second;% 1sec
%params.behav.UT_distance_from_CO=10; %m
params.behav.distnace_UT_from_ball=6;%m
params.behav.min_diff_UT=30;% samples
params.behav.dis_before_after_ut=1; %m
params.behav.time_before_ut_window=3; %s

params.behav.bins_to_remove_from_edge_CO_hist=1;
params.behav.manual_min_dis_from_CO=100; %for manual correction check that the CO is close
%new flight criteria (distance and time between bsp samples)
% when shuffling, we shuffle each flight separately. This is a proxy for when a new flight begins:
% 1. consecutive samples that have a jump in space
% 2. consecutive samples that have a time gap
% can be done with closer attention...
params.behav.dis_criteria = [20 20];
params.behav.new_flight_time_criteria = 1e5;

% interpolation of co with missing data
params.behav.dis_between_bats_interpolate_thresh=1; %(m)
params.behav.min_dist_to_co_to_interp=5; %m
params.behav.min_hole_size_to_interp=2; %sec


%save
behav_params=params.behav;
param_file_name=fullfile(param_folder,'behav_params.mat');
save(param_file_name, '-struct', 'behav_params')

%% parameters for solo analysis
%parameters for Tuning curve calculation 
% -  DO I NEED IT??? -------------
tunnel_limits=[0 135];
params.solo.solo_X_min= tunnel_limits(1);
params.solo.solo_X_max= tunnel_limits(2);
params.solo.solo_X_n_bins = params.solo.solo_X_max * 2; %0.5 meter bin
params.solo.solo_X_bin_size = (params.solo.solo_X_max-params.solo.solo_X_min)/params.solo.solo_X_n_bins;
params.solo.solo_X_bins_vector=params.solo.solo_X_min:params.solo.solo_X_bin_size:params.solo.solo_X_max;
params.solo.solo_X_bins_vector_of_centers=params.solo.solo_X_bins_vector(1:end-1)+params.solo.solo_X_bin_size/2;
params.solo.solo_time_spent_minimum_for_1D_bins=0.75;
params.solo.frames_per_second=params.behav.frame_per_second;
% for 2d comparison:
params.solo.allo_X_bin_size_2D=3; %changded to 3 meter
params.solo.allo_X_bins_vector_2D=tunnel_limits(1):params.solo.allo_X_bin_size_2D:tunnel_limits(2);
params.solo.allo_X_bins_vector_of_centers_2D=params.solo.allo_X_bins_vector_2D(1:end-1)+params.solo.allo_X_bin_size_2D/2;
params.solo.ker_SD=1.5;
%new flight criteria (distance and time between bsp samples)
% when shuffling, we shuffle each flight separately. This is a proxy for when a new flight begins:
% 1. consecutive samples that have a jump in space
% 2. consecutive samples that have a time gap
% can be done with closer attention...
params.solo.dis_criteria = params.behav.dis_criteria;
params.solo.new_flight_time_criteria = params.behav.new_flight_time_criteria;

%save
solo_params=params.solo;
param_file_name=fullfile(param_folder,'solo_params.mat');
save(param_file_name, '-struct', 'solo_params')

%% field detection params
params.fields.bin_centers=params.solo.solo_X_bins_vector_of_centers;
params.fields.bin_edges=params.solo.solo_X_bins_vector;
params.fields.bin_size = params.solo.solo_X_bin_size;
params.fields.bin_limits = [params.solo.solo_X_min params.solo.solo_X_max];
params.fields.ker_SD = params.solo.ker_SD;
params.fields.min_time_spent_per_bin = params.solo.solo_time_spent_minimum_for_1D_bins;
params.fields.ker_type = 'gaussian';
params.fields.shuffles_num = 1000;
params.fields.shuffles_max_shift = 30;

params.fields.FR_thr = 0.5;
params.fields.overlap_href = 0.5; % href=horizontal reference
params.fields.width_href = 0.2;
params.fields.width_prc = [5 95]; % TODO: consider chaning to 2.5-97.5 %
params.fields.min_spikes = 10;
params.fields.min_flights_with_spikes = 5;
params.fields.min_flights_with_spikes_prc = 0.2;
params.fields.local_shuffle.margin = 0.5; % relative to field width
params.fields.local_shuffle.n_shuffles = 1000;
params.fields.local_shuffle.max_shift = 30;
params.fields.local_shuffle.signif_SI_prc = 95;
params.fields.valid_speed_pos = [10 187.5];
params.fields.parmaset=1;
%save
fields_params=params.fields;
param_file_name=fullfile(param_folder,'fields_params.mat');
save(param_file_name, '-struct', 'fields_params')
%% solo shuffle:
params.solo_shuffle_params.shuffles_num=1000;
params.solo_shuffle_params.shuffles_max_shift=30; %check why with Tamir


params.solo_shuffle_params.cell_co_solo_initial_analysis_struct_folder = params.dirs.cell_co_solo_initial_analysis_struct_folder;
params.solo_shuffle_params.solo_shuffle_folder_name=params.dirs.solo_shuffle_folder_name;
%save
fields_params=params.solo_shuffle_params;
param_file_name=fullfile(param_folder,'solo_shuffle_params.mat');
save(param_file_name, '-struct', 'fields_params')

%% parameters for CO analysis basic analysis
params.co.combine_rectangle_min_sep=3;
params.co.min_coverage_of_ego_bin_2D=0.8;
params.co.time_spent_minimum_for_1D_bins=0.2;
params.co.frames_per_second=params.behav.frame_per_second;
params.co.alpha_val=5;
params.co.time_before_after_co_for_co_window=params.behav.time_before_after_co_for_co_window; %(seconds) for defining the window around the co
params.co.signif_by_SI=1; %flag on signif ego cell by zscore of SI
params.co.min_n_spike_ego_corr=20;
params.co.solo_X_bins_vector=params.solo.solo_X_bins_vector;

% a. Egocentric time bins: 
params.co.ker_SD=params.fields.ker_SD;
params.co.time_before_after_co=params.behav.time_before_after_co;
params.co.time_X_min= -params.co.time_before_after_co * 1e6;
params.co.time_X_max=params.co.time_before_after_co * 1e6;
params.co.time_n_bins = 40;  
params.co.time_X_bin_size=(params.co.time_X_max - params.co.time_X_min)/params.co.time_n_bins;
params.co.time_X_bins_vector=params.co.time_X_min:params.co.time_X_bin_size:params.co.time_X_max;
params.co.time_X_bins_vector_of_centers=params.co.time_X_bins_vector(1:end-1)+params.co.time_X_bin_size/2;

% b. Egocentric distance bins
params.co.dis_before_after_co=params.behav.dis_before_after_co;
params.co.dis_X_min = -params.co.dis_before_after_co;
params.co.dis_X_max = params.co.dis_before_after_co;
params.co.dis_n_bins = params.co.time_n_bins; %use the same number of bins like in time
params.co.dis_X_bin_size = (params.co.dis_X_max-params.co.dis_X_min)/params.co.dis_n_bins;
params.co.dis_X_bins_vector=params.co.dis_X_min:params.co.dis_X_bin_size:params.co.dis_X_max;
params.co.dis_X_bins_vector_of_centers=params.co.dis_X_bins_vector(1:end-1)+params.co.dis_X_bin_size/2;

%d. significant differnce between co and solo allocentric representation
params.co.sig_bins_width = 3;
params.co.min_flights_per_bin = 3;
params.co.tunnel_end=tunnel_limits(2);
%e. 2D tuning distance vs allocentric (trying few different bin sizes)
params.co.allo_X_bin_size_2D=params.solo.allo_X_bin_size_2D; %changded to 3 meter
params.co.dis_X_bin_size_2D=3; %changded to 3 meter
n_bins=round((params.co.dis_X_max-params.co.dis_X_min)/params.co.dis_X_bin_size_2D);
params.co.allo_X_bins_vector_2D=tunnel_limits(1):params.co.allo_X_bin_size_2D:tunnel_limits(2);
params.co.allo_X_bins_vector_of_centers_2D=params.co.allo_X_bins_vector_2D(1:end-1)+params.co.allo_X_bin_size_2D/2;
params.co.dis_X_bins_vector_2D=linspace(params.co.dis_X_min,params.co.dis_X_max,n_bins);
params.co.dis_X_bins_vector_of_centers_2D=params.co.dis_X_bins_vector_2D(1:end-1)+params.co.dis_X_bin_size_2D/2;
params.co.time_spent_minimum_for_2D_bins=0.2; 
params.co.sigma_a=1.5;
params.co.hsize=5*round(params.co.sigma_a)+1;
params.co.legalize_by_neighbor_bins_flag=1;
% for computing reduction in velocity:
params.co.long_dis_thresh=20;
params.co.short_dis_thresh=5;
%save
co_params=params.co;
param_file_name=fullfile(param_folder,'co_params.mat');
save(param_file_name, '-struct', 'co_params')
%% Per field params:
params.per_field.per_field_to_plot=1; %1= width based on half hight 2=width based on 5-95% 3=width based on all spikes
params.per_field.time_spent_minimum_for_1D_bins_per_field=0.1;
params.per_field.frames_per_second=100;
%bins params:
params.per_field.time_before_after_co=params.behav.time_before_after_co;
params.per_field.time_X_min= -params.per_field.time_before_after_co * 1e6;
params.per_field.time_X_max=params.per_field.time_before_after_co * 1e6;
params.per_field.time_per_field_n_bins = 30;  
params.per_field.time_per_field_X_bin_size=(params.per_field.time_X_max - params.per_field.time_X_min)/params.per_field.time_per_field_n_bins;
params.per_field.time_per_field_X_bins_vector=params.per_field.time_X_min:params.per_field.time_per_field_X_bin_size:params.per_field.time_X_max;
params.per_field.time_per_field_bin_vec_of_center=params.per_field.time_per_field_X_bins_vector(1:end-1)+params.per_field.time_per_field_X_bin_size/2;

% b. Egocentric distance bins
params.per_field.dis_before_after_co=params.behav.dis_before_after_co;
params.per_field.dis_X_min = -params.per_field.dis_before_after_co;
params.per_field.dis_X_max = params.per_field.dis_before_after_co;
params.per_field.dis_per_field_n_bins = params.per_field.time_per_field_n_bins; %use the same number of bins like in time
params.per_field.dis_per_field_X_bin_size = (params.per_field.dis_X_max-params.per_field.dis_X_min)/params.per_field.dis_per_field_n_bins;
params.per_field.dis_per_field_X_bins_vector=params.per_field.dis_X_min:params.per_field.dis_per_field_X_bin_size:params.per_field.dis_X_max;
params.per_field.dis_per_field_bin_vec_of_center=params.per_field.dis_per_field_X_bins_vector(1:end-1)+params.per_field.dis_per_field_X_bin_size/2;

params.per_field.n_bins=params.per_field.time_per_field_X_bins_vector;

%shuffle:
params.per_field.alpha_val=5;
params.per_field.num_shuffles_per_field=10000;
params.per_field.min_offset_perc=0.1; % perc of length of data need to think about this parameter!
params.per_field.benf_correct=10;

%width:
params.per_field.width_at_heigth=50;%width of per field
%smooth:
params.per_field.old_smooth=0;
params.per_field.smooth_window=3;
params.per_field.smooth_type=1;
params.per_field.smooth_tol=1;

%Per field population params:
%-----------------------------
%1. duing all pop analysis on nice place cells: for now as in pop vec analysis
params.per_field.min_n_spike=100;  
params.per_field.SI_threshold=1;
params.per_field.buffer_for_solo_shuffle_per_field=2;%m
%2. duing all pop analysis on nice tuning per field:
params.per_field.min_n_spike_per_field=30;
params.per_field.min_r_length_per_field=params.co.min_coverage_of_ego_bin_2D; %min length of non nan firing rate in per field

%3. tuning width params:
params.per_field.min_dis_pos_neg=0.25*1e6;%us


%save
per_field_params=params.per_field;
param_file_name=fullfile(param_folder,'per_field_params.mat');
save(param_file_name, '-struct', 'per_field_params')

%% parameters for CO shuffle code
%firing rate parameters

% a. general parameters
params.co_shuffle.dis_criteria = params.behav.dis_before_after_co;
params.co_shuffle.time_spent_minimum_for_1D_bins=params.co.time_spent_minimum_for_1D_bins;
params.co_shuffle.frames_per_second=params.co.frames_per_second;
params.co_shuffle.alpha_val=5; %

% b. egocentric parameters
params.co_shuffle.dis_X_bins_vector_of_centers=params.co.dis_X_bins_vector_of_centers;
params.co_shuffle.dis_X_bins_vector=params.co.dis_X_bins_vector;

% c. coherence parameters (wider bins)
params.co_shuffle.coherence_X_min = params.co.dis_X_min;
params.co_shuffle.coherence_X_max = -params.co.dis_X_min;
params.co_shuffle.coherence_n_bins = 20;
params.co_shuffle.coherence_X_bin_size = (params.co_shuffle.coherence_X_max-params.co_shuffle.coherence_X_min)/params.co_shuffle.coherence_n_bins;
params.co_shuffle.coherence_X_bins_vector=params.co_shuffle.coherence_X_min:params.co_shuffle.coherence_X_bin_size:params.co_shuffle.coherence_X_max;
params.co_shuffle.coherence_X_bins_vector_of_centers=params.co_shuffle.coherence_X_bins_vector(1:end-1)+params.co_shuffle.coherence_X_bin_size/2;

% d. allocentric parameters
params.co_shuffle.bin_centers=params.solo.solo_X_bins_vector_of_centers;
params.co_shuffle.allo_bin_size = params.solo.solo_X_bin_size;
params.co_shuffle.allo_bin_limits = [params.solo.solo_X_min params.solo.solo_X_max];
params.co_shuffle.ker_SD = 1.5;
params.co_shuffle.min_time_spent_per_meter = params.solo.solo_time_spent_minimum_for_1D_bins;
params.co_shuffle.ker_type = 'gaussian';
params.co_shuffle.pos_X_bins_vector_of_centers=params.co_shuffle.bin_centers;


%shuffling parameters
params.co_shuffle.co_shuffle_folder_name = params.dirs.co_shuffle_folder_name;

params.co_shuffle.n_shuffles =1000;
params.co_shuffle.n_shuffles = params.co_shuffle.n_shuffles+1;
params.co_shuffle.ego_shuffle = 1; % shuffle egocentric distances, not only allocentric position
params.co_shuffle.cell_co_solo_initial_analysis_struct_folder = params.dirs.cell_co_solo_initial_analysis_struct_folder;

%save
co_shuffle_params=params.co_shuffle;
param_file_name=fullfile(param_folder,'co_shuffle_params.mat');
save(param_file_name, '-struct', 'co_shuffle_params')

%% inclusion:
params.inclusion_cells.max_for_pyramidal=5;
params.inclusion_cells.min_n_spike_in_air=100; % per dir
params.inclusion_cells.min_co_per_dir=10; %for now need to check it!
%signif place cell:
params.inclusion_cells.SI_thr_solo=1;
params.inclusion_cells.SI_thr_shuffle_solo=95; %Tamir had it on 99 decide what we want

% signif egocentric cell:
params.inclusion_cells.SI_thr_shuffle_co=params.inclusion_cells.SI_thr_shuffle_solo; %
params.inclusion_cells.min_spikes_during_co = 30; 
params.inclusion_cells.SI_thr_co = 0.1;
params.inclusion_cells.min_even_odd = 0.2;

%dirs:
params.inclusion_cells.solo_shuffle_folder_name=params.dirs.solo_shuffle_folder_name;
params.inclusion_cells.co_shuffle_folder_name=params.dirs.co_shuffle_folder_name;

params.inclusion_cells.cell_co_solo_initial_analysis_struct_folder = params.dirs.cell_co_solo_initial_analysis_struct_folder;
params.inclusion_cells.inclusion_cells_folder_name= params.dirs.inclusion_cells_folder_name;
%save
inclusion_cells_params=params.inclusion_cells;
param_file_name=fullfile(param_folder,'inclusion_cells_params.mat');
save(param_file_name, '-struct', 'inclusion_cells_params')


%% CO population params:

params.co_population.min_spikes = params.inclusion_cells.min_spikes_during_co; %during co per direction
params.co_population.min_ego_inf = params.inclusion_cells.SI_thr_co;
params.co_population.min_even_odd = params.inclusion_cells.min_even_odd;
params.co_population.alpha_val = params.inclusion_cells.SI_thr_shuffle_co;
params.co_population.interneuron_firing_rate = 10; %??????
params.co_population.max_for_pyramidal=params.inclusion_cells.max_for_pyramidal;
%save
co_population_params=params.co_population;
param_file_name=fullfile(param_folder,'co_population_params.mat');
save(param_file_name, '-struct', 'co_population_params')


%% population vector params:
params.population_vector.n_shuffles=100;%temp
% controls for behavior:
%data in y:
params.population_vector.run_only_on_intersent_data=0;
params.population_vector.vel_reduction_thresh=0.9;
params.population_vector.run_only_fast_co=1;
%cell selection:
params.population_vector.SI_threshold=params.per_field.SI_threshold;
params.population_vector.min_n_spike=params.per_field.min_n_spike;
params.population_vector.solo_time_spent_minimum_for_1D_bins=0.03;%need to check!

%egocentric bins:
params.population_vector.ego_range=[-40 40];
params.population_vector.window=10;
params.population_vector.step=2;
params.population_vector.ego_bin_start=params.population_vector.ego_range(1):params.population_vector.step:params.population_vector.ego_range(2)-params.population_vector.window;
params.population_vector.ego_bin_end=params.population_vector.ego_range(1)+params.population_vector.window:params.population_vector.step:params.population_vector.ego_range(2);
params.population_vector.ego_bin_center=[params.population_vector.ego_bin_start+params.population_vector.ego_bin_end]./2;
params.population_vector.ego_bin_dis=params.population_vector.ego_range(1):params.population_vector.step:params.population_vector.ego_range(2);

%params.population_vector.full_data=0;

%allocentric tuning params:
params.population_vector.allo_X_bin_size = 1;%meter
params.population_vector.allo_X_min= tunnel_limits(1);
params.population_vector.allo_X_max= tunnel_limits(2);
params.population_vector.allo_X_n_bins = params.solo.solo_X_max; %1 meter bin
params.population_vector.allo_X_bin_size = (params.population_vector.allo_X_max-params.population_vector.allo_X_min)/params.population_vector.allo_X_n_bins;
params.population_vector.allo_X_bins_vector=params.population_vector.allo_X_min:params.population_vector.allo_X_bin_size:params.population_vector.allo_X_max;
params.population_vector.allo_X_bins_vector_of_centers=params.population_vector.allo_X_bins_vector(1:end-1)+params.population_vector.allo_X_bin_size/2;
params.population_vector.allo_time_spent_minimum_for_1D_bins=0.03;
params.population_vector.frames_per_second=params.behav.frame_per_second;
params.population_vector.ker_SD=params.fields.ker_SD;


%save
population_vector_params=params.population_vector;
param_file_name=fullfile(param_folder,'population_vector_params.mat');
save(param_file_name, '-struct', 'population_vector_params')

