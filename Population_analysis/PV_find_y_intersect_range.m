
function   [range_intersect, intersect_by_union]=PV_find_y_intersect_range(data,bin_i,allo_X_bins_vector,dir_i)

solo_per_bin_bsp_ind=(data(dir_i).solo.bsp.flight_x_pos>=allo_X_bins_vector(bin_i)& data(dir_i).solo.bsp.flight_x_pos<allo_X_bins_vector(bin_i+1));
co_per_bin_bsp_ind=(data(dir_i).co.bsp.flight_x_pos>=allo_X_bins_vector(bin_i)& data(dir_i).co.bsp.flight_x_pos<allo_X_bins_vector(bin_i+1));
if sum(solo_per_bin_bsp_ind(:))~=0 & sum(co_per_bin_bsp_ind(:))~=0
    
    % get y data in bin:
    %--------------------
    solo_y_bsp_pos_per_x_bin=data(dir_i).solo.bsp.flight_y_pos(solo_per_bin_bsp_ind);
    co_y_bsp_pos_per_x_bin=data(dir_i).co.bsp.flight_y_pos(co_per_bin_bsp_ind);
    
    %compute intersection
    %--------------------
    range_intersect=[nanmax([min(solo_y_bsp_pos_per_x_bin),min(co_y_bsp_pos_per_x_bin)]),nanmin([max(solo_y_bsp_pos_per_x_bin),max(co_y_bsp_pos_per_x_bin)])];
    
    % compute intersect by union:
    %----------------------------
    union_range=[nanmin([min(solo_y_bsp_pos_per_x_bin),min(co_y_bsp_pos_per_x_bin)]),nanmax([max(solo_y_bsp_pos_per_x_bin),max(co_y_bsp_pos_per_x_bin)])];
    intersect_by_union=(range_intersect(2)-range_intersect(1))/(union_range(2)-union_range(1));
    
    %check if there is no intersection
    %---------------------------------
    if min(solo_y_bsp_pos_per_x_bin)>max(co_y_bsp_pos_per_x_bin) | min(co_y_bsp_pos_per_x_bin)>max(solo_y_bsp_pos_per_x_bin)
        intersect_by_union=0;
        range_intersect=[nan nan];
    end
    
    
else
    range_intersect=[nan nan];
    intersect_by_union=nan;
end