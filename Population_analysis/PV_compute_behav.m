function   [intersect_by_union,y_dev_per_pos,vel]=PV_compute_behav(data,allo_X_bins_vector,ego_bin_start,ego_bin_end,dir_i)
 %take relevent data of co:
co_y=data(dir_i).co.bsp.flight_y_pos';
co_y=co_y(:);co_y(isnan(co_y))=[];
co_deviation_from_solo=nan*zeros(length(co_y),1);
co_x=data(dir_i).co.bsp.flight_x_pos';
co_x=co_x(:);co_x(isnan(co_x))=[];
co_ts=data(dir_i).co.bsp.flight_ts';
co_ts=co_ts(:);co_ts(isnan(co_ts))=[];
dis_x=data(dir_i).co.bsp.flight_dis';
dis_x=dis_x(:);dis_x(isnan(dis_x))=[];

%intialize mat:
co_deviation_from_solo=nan*zeros(length(co_y),1);
intersect_by_union=nan*zeros(length(allo_X_bins_vector)-1,1);
y_dev_per_pos=nan*zeros(length(ego_bin_start),1);
vel=nan*zeros(length(ego_bin_start),1);

%% go over allo bins to compute y deviation and intersect/union:


for allo_bin_i=1:length(allo_X_bins_vector)-1
    %find solo relevant for allo bin:
    solo_allo_bin_ind=(data(dir_i).solo.bsp.flight_x_pos>=allo_X_bins_vector(allo_bin_i)& data(dir_i).solo.bsp.flight_x_pos<allo_X_bins_vector(allo_bin_i+1));
    solo_y=data(dir_i).solo.bsp.flight_y_pos;
    solo_y_per_allo_bin=solo_y(solo_allo_bin_ind);
    if ~isempty(solo_y_per_allo_bin)
        mean_solo_y_per_x_pos=nanmean(solo_y_per_allo_bin);
        % find co relevant for allo bin:
        co_allo_bin_ind=(co_x>=allo_X_bins_vector(allo_bin_i)& co_x<allo_X_bins_vector(allo_bin_i+1));
        %compute deviation of y in co vs mean solo per pos:
        co_deviation_from_solo(co_allo_bin_ind)=abs(co_y(co_allo_bin_ind)'-mean_solo_y_per_x_pos);
        % compute intersection by union:
        [~, intersect_by_union(allo_bin_i)]=PV_find_y_intersect_range(data,allo_bin_i,allo_X_bins_vector,dir_i);
        
    end
end
%% go over egocentric bins to compute velocity and dev in y per ago bin:
for ego_bin_i=1:length(ego_bin_start)
    % find relevant inds for solo and co:
    co_ego_bin_ind=(dis_x>=ego_bin_start(ego_bin_i)& dis_x<ego_bin_end(ego_bin_i));
    if sum(co_ego_bin_ind(:))~=0
        %1. compute mean velocity:
        factor_to_sec=1e6;
        x=co_x(co_ego_bin_ind);
        y=co_y(co_ego_bin_ind);
        ts=co_ts(co_ego_bin_ind);
        vel(ego_bin_i)=nanmean(compute_velocity(x,y,ts,factor_to_sec));
        
        %2. compute mean deviation in y per ego bin:
        y_dev_per_pos(ego_bin_i)=nanmean(co_deviation_from_solo(co_ego_bin_ind));

    else
        vel(ego_bin_i)=nan;
        y_dev_per_pos(ego_bin_i)=nan;
        
    end
end