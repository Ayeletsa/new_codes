%temp code to plot co position distance from ball position:
clear
dir_day_struct='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\go_over_behavior_again\';
dir_infoday_struct=dir(dir_day_struct);
ball_pos_dir='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\ball_pos\';
dir_info_ball_pos=dir(ball_pos_dir);
for day_i=3:length(dir_infoday_struct)
    load(fullfile(dir_day_struct,dir_infoday_struct(day_i).name));
    %load ball pos
    day=dir_infoday_struct(day_i).name(end-11:end-4);
    
    ind_day_ball=find(~cellfun(@isempty,regexp({dir_info_ball_pos.name},day)));
    ind_day_ball=ind_day_ball(1);
    load(fullfile(ball_pos_dir,dir_info_ball_pos(ind_day_ball).name));
    
    
    behav_file_name=dir_infoday_struct(day_i).name;
     bat=behav_file_name(22:25);
    general_behavior_data_file_name=['D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\general_behave_analysis_bat',bat,'_day',day,'.mat\'];
    load(general_behavior_data_file_name)
    
    % find co pos:
    co_ind=behavioral_modes.CO_point ;
    bat_pos_during_co=bsp_proc_data(tag_i).pos(co_ind,1);  
    % find distance between ball and co
    dist_to_balls=abs([bat_pos_during_co-ball_1_pos(1),bat_pos_during_co-ball_2_pos(1)]);
    % find the minimum of the two distances
    min_dist_to_ball=min(dist_to_balls,[],2);
    %compute hist of ball pos
    figure
    x_vec=1:2:71;
    hist(min_dist_to_ball,x_vec)
    n_co_lower_than_10=length(find(min_dist_to_ball<=10));
    title(sprintf('distance from ball - bat %s day %s\n #co lower than 10 = %d',bat,day,n_co_lower_than_10))
    %save figures:
    fig_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\co_dist_from_ball\';
    fig_name=fullfile(fig_folder,sprintf('distance from ball - bat %s day %s.jpg',bat,day));
    saveas(gcf,fig_name)
    
    close all
end