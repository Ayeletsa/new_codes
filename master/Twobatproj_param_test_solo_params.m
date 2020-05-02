function Twobatproj_param_test_solo_params(param_folder,kernel,cells_struct_dir,FR_thr)
%for testig- REMOVE!!
%kernel=1.5;
%% folders:
% data input:
%params.dirs.cells_struct_dir='L:\Data\2batproj\Data_Nlg_Proc\yr_2019_bat_2336\cell_structs\';
%params.dirs.cells_struct_dir='D:\Ayelet\Data\Data_Nlg_Proc\yr_2018_bat_2389\cell_structs\';
params.dirs.cells_struct_dir=cells_struct_dir;
% data output:
main_analysis_dir='D:\Ayelet\2bat_proj\Analysis\new_code\';
params.dirs.behave_day_struct_folder=[main_analysis_dir,'\analysis_structs\behavioral_modes\day_structs\'];
params.dirs.ball_position_folder=[main_analysis_dir,'\analysis_structs\behavioral_modes\ball_pos\'];
%params.dirs.cell_co_solo_initial_analysis_struct_folder=[main_analysis_dir,'\analysis_structs\co_solo_initial_analysis\'];
params.dirs.cell_co_solo_initial_analysis_struct_folder=[main_analysis_dir,'\analysis_structs\co_solo_initial_analysis_k_',num2str(kernel),'_th_',num2str(FR_thr),'\'];
if ~exist(params.dirs.cell_co_solo_initial_analysis_struct_folder, 'dir')
   mkdir(params.dirs.cell_co_solo_initial_analysis_struct_folder) 
end
params.dirs.co_shuffle_folder_name = [main_analysis_dir,'\analysis_structs\co_shuffling_struct'];
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
params.behav.min_solo_length=2*params.behav.frame_per_second; %samples
params.behav.min_tracking_length=3*params.behav.frame_per_second; %samples
params.behav.dist_thresh_tracking=20; %meters
params.behav.dist_thresh_solo=25; %meters
params.behav.dist_thresh_CO=5; %meters
params.behav.CO_window=[20 20]; %meters
params.behav.max_wind_CO=2*params.behav.frame_per_second;
params.behav.min_time_before_CO=2*params.behav.frame_per_second;
params.behav.UT_window=params.behav.frame_per_second;
params.behav.UT_time_from_CO=2*params.behav.frame_per_second;% 1sec
params.behav.UT_distance_from_CO=10; %m
params.behav.bins_to_remove_from_edge_CO_hist=1;
params.behav.manual_min_dis_from_CO=100; %for manual correction check that the CO is close
params.behav.frame_per_second=100;


%save
behav_params=params.behav;
param_file_name=fullfile(param_folder,'behav_params.mat');
save(param_file_name, '-struct', 'behav_params')

%% parameters for solo analysis
%parameters for Tuning curve calculation

params.solo.solo_X_min= 0;
params.solo.solo_X_max= 135;
params.solo.solo_X_n_bins = params.solo.solo_X_max * 2;
params.solo.solo_X_bin_size = (params.solo.solo_X_max-params.solo.solo_X_min)/params.solo.solo_X_n_bins;
params.solo.solo_X_bins_vector=params.solo.solo_X_min:params.solo.solo_X_bin_size:params.solo.solo_X_max;
params.solo.solo_X_bins_vector_of_centers=params.solo.solo_X_bins_vector(1:end-1)+params.solo.solo_X_bin_size/2;
params.solo.solo_time_spent_minimum_for_1D_bins=0.75;
params.solo.frames_per_second=100;
%
% field_detection_bin_size = 0.5; %m
% field_detection_X_bins_vector=solo_X_min:field_detection_bin_size:solo_X_max;
% field_detection_X_bins_vector_of_centers=field_detection_X_bins_vector(1:end-1)+field_detection_bin_size/2;
%
params.solo.field_detection_X_bins_vector=params.solo.solo_X_bins_vector;

%new flight criteria (distance and time between bsp samples)
% when shuffling, we shuffle each flight separately. This is a proxy for when a new flight begins:
% 1. consecutive samples that have a jump in space
% 2. consecutive samples that have a time gap
% can be done with closer attention...
params.solo.dis_criteria = [20 20];
params.solo.new_flight_time_criteria = 1e5;

% params for field detection code (find_fields_PSTH_Basic_withShuffle)
params.solo.ref_height_for_overlap = 0.5;
params.solo.minPeakFR = 0.5; %Hz
params.solo.min_spikes_perField = 10;
params.solo.percentilesTH = 85; %the firing rate TH has to be higher than 85% of the bins.
params.solo.stabilityTH = 0.15; % cell has to fire within the field (defind by the
% ref_height_for_width) in at least 15% of flights.
params.solo.num_of_flights_minimum = 5; % cell has to fire within the field (defind by the
% ref_height_for_width) in at least 5 flights.
params.solo.valley_2_fields_unite_TH = 1; %if the valley between two close fields does
% not go below 1Hz - unite them to one field.
params.solo.bin_edges=params.solo.field_detection_X_bins_vector;
params.solo.ref_height_for_width=0.2;
params.solo.ref_height_for_overlap=0.5;

%save
solo_params=params.solo;
param_file_name=fullfile(param_folder,'solo_params.mat');
save(param_file_name, '-struct', 'solo_params')

%% field detection params
params.fields.bin_size = params.solo.solo_X_bin_size;
params.fields.bin_limits = [params.solo.solo_X_min params.solo.solo_X_max];
params.fields.ker_SD = kernel;
params.fields.min_time_spent_per_meter = params.solo.solo_time_spent_minimum_for_1D_bins;
params.fields.ker_type = 'gaussian';
params.fields.shuffles_num = 1000;
params.fields.shuffles_max_shift = 30;

params.fields.FR_thr = FR_thr;
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

%% parameters for CO analysis basic analysis
params.co.time_spent_minimum_for_1D_bins=0.2;
params.co.time_spent_minimum_for_1D_bins_per_field=0.1;
params.co.frames_per_second=100;
params.co.num_shuffles_per_field=1000;
params.co.alpha_val=5;
params.co.manual_min_dis_from_CO=100; %for manual correction check that the CO is close
% a. relative time to co
params.co.time_before_after_co=params.behav.time_before_after_co;
params.co.dis_before_after_co=params.behav.dis_before_after_co;
params.co.time_X_min= -params.co.time_before_after_co * 1e6;
params.co.time_X_max=params.co.time_before_after_co * 1e6;
params.co.time_n_bins = 30;
params.co.time_X_bin_size=(params.co.time_X_max - params.co.time_X_min)/params.co.time_n_bins;
params.co.time_X_bins_vector=params.co.time_X_min:params.co.time_X_bin_size:params.co.time_X_max;
params.co.time_X_bins_vector_of_centers=params.co.time_X_bins_vector(1:end-1)+params.co.time_X_bin_size/2;

% b. distance between bats
params.co.dis_X_min = -params.co.dis_before_after_co;
params.co.dis_X_max = params.co.dis_before_after_co;
params.co.dis_n_bins = params.co.time_n_bins;
params.co.dis_X_bin_size = (params.co.dis_X_max-params.co.dis_X_min)/params.co.dis_n_bins;
params.co.dis_X_bins_vector=params.co.dis_X_min:params.co.dis_X_bin_size:params.co.dis_X_max;
params.co.dis_X_bins_vector_of_centers=params.co.dis_X_bins_vector(1:end-1)+params.co.dis_X_bin_size/2;

%c. allocentric during CO
params.co.allo_X_min= 0;
params.co.allo_X_max= 135;
params.co.allo_X_n_bins = params.co.allo_X_max * 2;
params.co.allo_X_bin_size = (params.co.allo_X_max-params.co.allo_X_min)/params.co.allo_X_n_bins;
params.co.allo_X_bins_vector=params.co.allo_X_min:params.co.allo_X_bin_size:params.co.allo_X_max;
params.co.allo_X_bins_vector_of_centers=params.co.allo_X_bins_vector(1:end-1)+params.co.allo_X_bin_size/2;

%d. significant differnce between co and solo allocentric representation
params.co.sig_bins_width = 3;
params.co.min_flights_per_bin = 3;



%e. 2D tuning distance vs allocentric (trying few different bin sizes)
params.co.allo_X_bin_size_2D=3; %changded to 3 meter
params.co.dis_X_bin_size_2D=3; %changded to 3 meter
params.co.time_spent_minimum_for_2D_bins=0.2; %test this
params.co.sigma_a=1.5;
params.co.hsize=5*round(params.co.sigma_a)+1;
params.co.legalize_by_neighbor_bins_flag=1;

% params
params.co.n_time_spent_bins = 8;
params.co. tunnel_end = 135;
params.co.time_spent_criteria = params.co.time_spent_minimum_for_1D_bins/2;

% ego and allo bins edges
params.co.ego_bins_width = params.co.dis_before_after_co*2/params.co.n_time_spent_bins;
params.co.ego_bins_edges = -params.co.dis_before_after_co:params.co.ego_bins_width:params.co.dis_before_after_co;
params.co.allo_bins_width = params.co.tunnel_end/params.co.n_time_spent_bins;
params.co.allo_bins_edges = params.co.tunnel_end:-params.co.allo_bins_width:0;
%per field
params.co.width_at_heigth=50;%width of per field
%save
co_params=params.co;
param_file_name=fullfile(param_folder,'co_params.mat');
save(param_file_name, '-struct', 'co_params')


%% parameters for CO shuffle code
%firing rate parameters

% a. general parameters
params.co_shuffle.dis_criteria = params.behav.dis_before_after_co;
params.co_shuffle.time_spent_minimum_for_1D_bins=params.co.time_spent_minimum_for_1D_bins;
params.co_shuffle.frames_per_second=params.co.frames_per_second;
params.co_shuffle.alpha_val=5; %

% b. egocentric parameters
params.co_shuffle.dis_X_min = params.co.dis_X_min;
params.co_shuffle.dis_X_max = -params.co.dis_X_min;
params.co_shuffle.dis_n_bins = 40;
params.co_shuffle.dis_X_bin_size = (params.co_shuffle.dis_X_max-params.co_shuffle.dis_X_min)/params.co_shuffle.dis_n_bins;
params.co_shuffle.dis_X_bins_vector=params.co_shuffle.dis_X_min:params.co_shuffle.dis_X_bin_size:params.co_shuffle.dis_X_max;
params.co_shuffle.dis_X_bins_vector_of_centers=params.co_shuffle.dis_X_bins_vector(1:end-1)+params.co_shuffle.dis_X_bin_size/2;

% c. coherence parameters (wider bins)
params.co_shuffle.coherence_X_min = params.co.dis_X_min;
params.co_shuffle.coherence_X_max = -params.co.dis_X_min;
params.co_shuffle.coherence_n_bins = 20;
params.co_shuffle.coherence_X_bin_size = (params.co_shuffle.coherence_X_max-params.co_shuffle.coherence_X_min)/params.co_shuffle.coherence_n_bins;
params.co_shuffle.coherence_X_bins_vector=params.co_shuffle.coherence_X_min:params.co_shuffle.coherence_X_bin_size:params.co_shuffle.coherence_X_max;
params.co_shuffle.coherence_X_bins_vector_of_centers=params.co_shuffle.coherence_X_bins_vector(1:end-1)+params.co_shuffle.coherence_X_bin_size/2;

% d. allocentric parameters
params.co_shuffle.tunnel_limits = [0 135];
params.co_shuffle.pos_X_min = params.co_shuffle.tunnel_limits(1);
params.co_shuffle.pos_X_max = params.co_shuffle.tunnel_limits(2);
params.co_shuffle.pos_n_bins = 240;
params.co_shuffle.pos_X_bin_size = (params.co_shuffle.pos_X_max-params.co_shuffle.pos_X_min)/params.co_shuffle.pos_n_bins;
params.co_shuffle.pos_X_bins_vector=params.co_shuffle.pos_X_min:params.co_shuffle.pos_X_bin_size:params.co_shuffle.pos_X_max;
params.co_shuffle.pos_X_bins_vector_of_centers=params.co_shuffle.pos_X_bins_vector(1:end-1)+params.co_shuffle.pos_X_bin_size/2;


%% shuffling parameters
params.co_shuffle.co_shuffle_folder_name = params.dirs.co_shuffle_folder_name;

params.co_shuffle.n_shuffles =1000;
params.co_shuffle.n_shuffles = params.co_shuffle.n_shuffles+1;
params.co_shuffle.ego_shuffle = 1; % shuffle egocentric distances, not only allocentric position
params.co_shuffle.cell_co_solo_initial_analysis_struct_folder = params.dirs.cell_co_solo_initial_analysis_struct_folder;

%save
co_shuffle_params=params.co_shuffle;
param_file_name=fullfile(param_folder,'co_shuffle_params.mat');
save(param_file_name, '-struct', 'co_shuffle_params')


%% CO population params:

params.co_population.min_spikes = 30;
params.co_population.min_ego_inf = 0.1;
params.co_population.min_even_odd = 0.2;
params.co_population.alpha = .05;
params.co_population.interneuron_firing_rate = 10;
params.co_population.max_for_pyramidal=5;
%save
co_population_params=params.co_population;
param_file_name=fullfile(param_folder,'co_population_params.mat');
save(param_file_name, '-struct', 'co_population_params')



