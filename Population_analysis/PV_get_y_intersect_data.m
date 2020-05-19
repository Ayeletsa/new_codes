
function data=PV_get_y_intersect_data(data,allo_X_bins_vector,allo_X_bins_vector_of_centers)
% run over y bins and for each bin look for y range that is intersect
% between solo and co get data only within this range:
for dir_i=1:2
    ind_intersect=struct();
    other_dir_ind_intersect=struct();
    
    for bin_i=1:length(allo_X_bins_vector_of_centers)
        range_intersect=PV_find_y_intersect_range(data,bin_i,allo_X_bins_vector,dir_i);
        
        %if sum(isnan(range_intersect))~=2
            % find ind out of intersection:
            %---------------------------------
            ind_intersect=PV_find_ind_data_out_of_intersect_y(ind_intersect,data,range_intersect,allo_X_bins_vector,bin_i,dir_i);
%         else
%             % if there is no intersection 
       % end
    end
    %% remove data:
    data=PV_remove_non_intersect_data(data,dir_i,ind_intersect);
    
end




