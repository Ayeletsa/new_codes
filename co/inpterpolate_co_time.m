function   co_time_usec=inpterpolate_co_time(bsp_proc_data,tag_i,co_ind,co_time_usec,dis_between_bats_interpolate_thresh,bsp_x_pos,other_bat_x_pos,bsp_ts,bsp_dis_m,time_before_after_co_for_co_window,dis_before_after_co,us_factor,frames_per_second,bat,day,ii_co,ii_dir,min_dist_to_co_to_interp,min_hole_size_to_interp)

co_ts_before=co_time_usec;
distance_between_bats_in_co=bsp_proc_data(tag_i).pos(co_ind,1)' - bsp_proc_data(3-tag_i).pos(co_ind,1)';
%check if the distance between bats is too large in co point:
if abs(distance_between_bats_in_co)>dis_between_bats_interpolate_thresh
    %1. co data:
    bsp_relative_times = bsp_ts - co_time_usec;
    bsp_time_criteria_ind = ( - time_before_after_co_for_co_window*us_factor < bsp_relative_times & bsp_relative_times < time_before_after_co_for_co_window*us_factor);
    bsp_dis_criteria_ind = ( abs(bsp_dis_m)<dis_before_after_co)';
    bsp_ind = bsp_time_criteria_ind & bsp_dis_criteria_ind;
    bsp_x_pos_at_co_to_interp = bsp_x_pos(bsp_ind);
    bsp_ts_usec_co_with_holes =  bsp_ts(bsp_ind)';
    other_bat_x_pos_to_interp =  other_bat_x_pos(bsp_ind)';
    dis_between_bats=bsp_x_pos_at_co_to_interp-other_bat_x_pos_to_interp;
    %2. find the hole:
    usual_ts_diff=us_factor/frames_per_second;
    ts_diff=diff(bsp_ts_usec_co_with_holes);
    ts_hole_logical=ts_diff>usual_ts_diff;
    ts_hole_start_ind=find(ts_hole_logical);
    ts_hole_end_ind=ts_hole_start_ind+1;
    if ~isempty(ts_hole_start_ind)
        ts_to_interpolate=[bsp_ts_usec_co_with_holes(1:ts_hole_start_ind(1))];
        for hole_i=1:length(ts_hole_start_ind)
            ts_to_interpolate=[ts_to_interpolate,bsp_ts_usec_co_with_holes(ts_hole_start_ind(hole_i)):usual_ts_diff:bsp_ts_usec_co_with_holes(ts_hole_end_ind(hole_i))];
            
        end
        ts_to_interpolate=[ts_to_interpolate,bsp_ts_usec_co_with_holes(ts_hole_end_ind(end):end)];
    
    
    % intrpolate:
    %1. for self pos:
    interp_self_pos=interp1(bsp_ts_usec_co_with_holes,bsp_x_pos_at_co_to_interp,ts_to_interpolate);
    %2. for other pos:
    interp_other_pos=interp1(bsp_ts_usec_co_with_holes,other_bat_x_pos_to_interp,ts_to_interpolate);
    
    % find co point:
    dis_between_bats_interp=interp_self_pos-interp_other_pos;
    
    distnace_other_from_self_shifted_plus=[0, dis_between_bats_interp];
    distnace_other_from_self_shifted_minus=[dis_between_bats_interp, 0];
    
    distance_change_sign_interp=find(sign(distnace_other_from_self_shifted_minus.*distnace_other_from_self_shifted_plus)==-1)-1;
    
    %deinfe new co ts:
    new_co_ts=ts_to_interpolate(distance_change_sign_interp);
    new_co_pos_self=interp_self_pos(distance_change_sign_interp);
    % find hole next to co:
    hole_ts=[bsp_ts_usec_co_with_holes(ts_hole_start_ind);bsp_ts_usec_co_with_holes(ts_hole_end_ind)];
    hole_ind=find(hole_ts(1,:)<new_co_ts & hole_ts(2,:)>new_co_ts);
    hole_ts_size=(hole_ts(2,hole_ind)-hole_ts(1,hole_ind))/us_factor; %sec
    distance_before=abs(dis_between_bats(ts_hole_start_ind(hole_ind)));
    distance_after=abs(dis_between_bats(ts_hole_end_ind(hole_ind)));

    %plot
    figure;
    plot(ts_to_interpolate,interp_self_pos,'g.')
    hold on;
    plot(ts_to_interpolate,interp_other_pos,'r.')
    plot(new_co_ts,new_co_pos_self,'k*')
    plot([co_ts_before co_ts_before],[0 135],'k')
    plot(bsp_ts_usec_co_with_holes,bsp_x_pos_at_co_to_interp,'k.')
    plot(bsp_ts_usec_co_with_holes,other_bat_x_pos_to_interp,'k.')

   
    
    %% correct ts only if:
    % if hole size lower than:
    if hole_ts_size< min_hole_size_to_interp
        co_time_usec=new_co_ts;
          title(sprintf('correct:hole size %.1f sec, distance before %.1f m, after %.1f m',hole_ts_size,distance_before,distance_after))
    elseif  distance_before<min_dist_to_co_to_interp | distance_after<min_dist_to_co_to_interp
         title(sprintf('correct:hole size %.1f sec, distance before %.1f m, after %.1f m',hole_ts_size,distance_before,distance_after))
        co_time_usec=new_co_ts;
    else
       title(sprintf('DO not correct:hole size %.1f sec, distance before %.1f m, after %.1f m',hole_ts_size,distance_before,distance_after))
        co_time_usec=[];
    end
    
    %%
     fig_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\co_time_interp\';
    saveas(gcf,fullfile(fig_folder,sprintf('interp_co_bat_%d_day_%s_dir_%d_co_%d.jpg',bat,day,ii_dir,ii_co)))
    close all
    else
        co_time_usec=[];
    end
end
