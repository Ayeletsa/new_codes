function per_field=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data)

load(per_field_params_file_name)

for field_i=1:size(field_edges,2)
    spikes_vec=spikes.dis_m(find(spikes.x_pos>field_edges(1,field_i) & spikes.x_pos<field_edges(2,field_i)));
    bsp_vec=bsp.dis_m(find(bsp.x_pos>field_edges(1,field_i) & bsp.x_pos<field_edges(2,field_i)));  per_field(field_i).number_of_spikes_per_field=length(spikes_vec);
    if ~isempty(spikes_vec)
        
        %per field for distance:
        %=============================================
        %1. compute tuning:
        [timespent_binned_new_smooth, ~, ~, tuning_dis_x_fr_per_field, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth_new_again ...
            (bsp_vec,spikes_vec,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0,old_smooth,smooth_window,smooth_type,smooth_tol);
        %2. run shuffle:
        dis_signif_field=shuffling_per_field_analysis_new_smooth(bsp_vec,spikes_vec,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth);
        %3. compute properties:
        per_field=per_field_prop(per_field,field_i,tuning_dis_x_fr_per_field,timespent_binned_new_smooth);
        %4. save data to struct:
        per_field(field_i).dis_signif_field=dis_signif_field;
        per_field(field_i).tuning_dis_x_fr_per_field=tuning_dis_x_fr_per_field;
        
        
        %per field for time:
        %=============================================
        spikes_vec=spikes.time_to_co(find(spikes.x_pos>field_edges(1,field_i) & spikes.x_pos<field_edges(2,field_i)));
        bsp_vec=bsp.time_to_co(find(bsp.x_pos>field_edges(1,field_i) & bsp.x_pos<field_edges(2,field_i)));
                 
        %1. compute tuning:
        [~, ~, ~, tuning_time_fr_per_field, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth_new_again ...
            (bsp_vec,spikes_vec,time_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0,old_smooth,smooth_window,smooth_type,smooth_tol);
        %2. run shuffle:
        time_signif_field=shuffling_per_field_analysis_new_smooth(bsp_vec,spikes_vec,time_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth);
        %3. save data to struct:        
        per_field(field_i).time_signif_field=dis_signif_field;
        per_field(field_i).tuning_time_fr_per_field=tuning_time_fr_per_field;
        
        %per field solo variability:
        %===================================
        per_field(field_i).solo_firing_rates_time_bins=per_field_solo_variability(field_edges,solo_data,time_X_bin_size,frames_per_second);
                
    end
end