function [per_field,per_field_tunings_corrs]=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data)

load(per_field_params_file_name)

for field_i=1:size(field_edges,2)
    
    %% shuffle:
    %1. arange data to shuffle:
    %spikes:
    spike_dis_mat=nan*zeros(size(spikes.x_pos));
    spike_time_mat=nan*zeros(size(spikes.x_pos));
    spike_ts_mat=nan*zeros(size(spikes.x_pos));
    
    [r,c]=find(spikes.x_pos>field_edges(1,field_i) & spikes.x_pos<field_edges(2,field_i));
    spike_dis_mat(sub2ind(size(spike_dis_mat),r,c))=spikes.dis_m(sub2ind(size(spike_dis_mat),r,c));
    spike_time_mat(sub2ind(size(spike_time_mat),r,c))=spikes.time_to_co(sub2ind(size(spike_time_mat),r,c));
    spike_ts_mat(sub2ind(size(spike_ts_mat),r,c))=spikes.ts_usec(sub2ind(size(spike_ts_mat),r,c));
    %bsp:
    bsp_dis_mat=nan*zeros(size(bsp.x_pos));
    bsp_time_mat=nan*zeros(size(bsp.x_pos));
    bsp_ts_mat=nan*zeros(size(bsp.x_pos));
    bsp_pos_x_mat=nan*zeros(size(bsp.x_pos));
    
    [r,c]=find(bsp.x_pos>field_edges(1,field_i) & bsp.x_pos<field_edges(2,field_i));
    bsp_dis_mat(sub2ind(size(bsp_dis_mat),r,c))=bsp.dis_m(sub2ind(size(bsp_dis_mat),r,c));
    bsp_time_mat(sub2ind(size(bsp_time_mat),r,c))=bsp.time_to_co(sub2ind(size(bsp_time_mat),r,c));
    bsp_ts_mat(sub2ind(size(bsp_ts_mat),r,c))=bsp.ts_usec(sub2ind(size(bsp_ts_mat),r,c));
    bsp_pos_x_mat(sub2ind(size(bsp_pos_x_mat),r,c))=bsp.x_pos(sub2ind(size(bsp_pos_x_mat),r,c));
    
    % find relevant co for field:
    co_in_field=find(sum(~isnan(bsp_dis_mat),2));
    % remove non-relevant co:
    spike_for_shffule.dis_m=spike_dis_mat(co_in_field,:);
    spike_for_shffule.time_to_co=spike_time_mat(co_in_field,:);
    spike_for_shffule.ts_usec=spike_ts_mat(co_in_field,:);
    
    bsp_for_shffule.dis_m=bsp_dis_mat(co_in_field,:);
    bsp_for_shffule.time_to_co=bsp_time_mat(co_in_field,:);
    bsp_for_shffule.ts_usec=bsp_ts_mat(co_in_field,:);
    
    %%
    
    spikes_vec=spike_dis_mat(:);spikes_vec=spikes_vec(~isnan(spikes_vec));
    bsp_vec=bsp_dis_mat(:);  bsp_vec=bsp_vec(~isnan(bsp_vec));
    per_field(field_i).number_of_spikes_per_field=length(spikes_vec);
    [~, ~, ~, r, ~,~] ...
        = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spikes_vec,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0);
    [ind_length,~,~]=find_length_of_consecutive_ind(find(~isnan(r)),length(r));
    
    
    if per_field(field_i).number_of_spikes_per_field>min_n_spike_per_field && max(ind_length)>=min_r_length_per_field*length(r)
        % 2. run shuffle:
        shuffles = rigid_shuffle_of_trials(bsp_for_shffule,spike_for_shffule,num_shuffles_per_field,min_offset_perc);
        
        %per field for distance:
        %=============================================
        %1. compute tuning:
        distance=1;
        [timespent_binned_new_smooth, ~, ~, tuning_dis_x_fr_per_field, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0);
        
        %         [timespent_binned_new_smooth, ~, ~, tuning_dis_x_fr_per_field, ~,~] ...
        %             = fn_compute_generic_1D_tuning_new_smooth_new_again ...
        %             (bsp_vec,spikes_vec,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0,old_smooth,smooth_window,smooth_type,smooth_tol);
        %         %2. run shuffle:
        dis_signif_field=shuffling_per_field_analysis_new_smooth(bsp_for_shffule,shuffles,dis_per_field_X_bins_vector,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth,min_dis_pos_neg,distance,benf_correct);
        %3. compute properties:
        per_field=per_field_prop(per_field,field_i,tuning_dis_x_fr_per_field,timespent_binned_new_smooth);
        %4. save data to struct:
        per_field(field_i).dis_signif_field=dis_signif_field;
        per_field(field_i).tuning_dis_x_fr_per_field=tuning_dis_x_fr_per_field;
         [per_field(field_i).solo_interp_time_tuning_shuffle, per_field(field_i).solo_interp_dis_tuning_shuffle]=create_solo_shuffle_per_field(solo_data,bsp_pos_x_mat,bsp_time_mat,bsp_dis_mat,time_per_field_bin_vec_of_center,dis_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second,buffer_for_solo_shuffle_per_field);
        
        %per field for time:
        %=============================================
        distance=0;
        spikes_vec=spikes.time_to_co(find(spikes.x_pos>field_edges(1,field_i) & spikes.x_pos<field_edges(2,field_i)));
        bsp_vec=bsp.time_to_co(find(bsp.x_pos>field_edges(1,field_i) & bsp.x_pos<field_edges(2,field_i)));
        
        %1. compute tuning:
        [~, ~, ~, tuning_time_fr_per_field, ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (bsp_vec,spikes_vec,time_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second, 0,0,0);
        
        %2. run shuffle:
        time_signif_field=shuffling_per_field_analysis_new_smooth(bsp_for_shffule,shuffles,time_per_field_X_bins_vector,time_per_field_bin_vec_of_center, time_spent_minimum_for_1D_bins_per_field, frames_per_second,num_shuffles_per_field,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth,min_dis_pos_neg,distance,benf_correct);
        %3. save data to struct:
        per_field(field_i).time_signif_field=time_signif_field;
        per_field(field_i).tuning_time_fr_per_field=tuning_time_fr_per_field;
        
        %per field solo variability:
        %===================================
        per_field(field_i).solo_firing_rates_time_bins=per_field_solo_variability(field_edges,solo_data,time_per_field_X_bin_size,frames_per_second);
    else
        per_field(field_i).tuning_time_fr_per_field=nan*zeros(1,length(time_per_field_bin_vec_of_center));
        per_field(field_i).tuning_dis_x_fr_per_field=nan*zeros(1,length(time_per_field_bin_vec_of_center));
    end
end
%% compute corrs per cell:
% run only if there is more than one field
if size(field_edges,2)>1
    per_field_tunings_corrs=corr(reshape([per_field.tuning_dis_x_fr_per_field],length(dis_per_field_bin_vec_of_center),[]),'rows','pairwise');
else
    per_field_tunings_corrs=nan;
end
