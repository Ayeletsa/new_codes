%% shuffling for per field analysis
%-------------------------------------------------------------------------
%here we are doing non rigid shuffling, meaning that the spikes are not in
%the same temporal relation in each shuffle.
% we just take all spikes and randomly distribute them between all
% positions

function [signif_field]=shuffling_per_field_analysis_new_smooth(bsp_vec,spikes_vec,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D, frames_per_second,num_shuffles,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth);
number_of_spikes=length(spikes_vec);

behav_length=length(bsp_vec);
parfor shuffle_i=1:num_shuffles
    if shuffle_i==1
        spike_vec_shuf=spikes_vec; %original data
    else
        %choosing random ind for the spikes
        rand_pos=randi(behav_length,number_of_spikes,1);
        %distribute the spikes new pos
        spike_vec_shuf=bsp_vec(rand_pos);
    end
    %compute tuning curve
    [~, ~, ~, r(shuffle_i,:), ~,~] ...
        = fn_compute_generic_1D_tuning_new_smooth_new_again ...
        (bsp_vec,spike_vec_shuf,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D, frames_per_second, 0,0,0,old_smooth,smooth_window,smooth_type,smooth_tol)
    
    std_r=nanstd(r(shuffle_i,:));
    mean_r=nanmean(r(shuffle_i,:));
    cv(shuffle_i)=std_r/mean_r;
    
    
    
end

%% check signif
%1.check if signif based on cv:
prct_alpha=prctile(cv,100-alpha_val);
if cv(1)>prct_alpha
    signif_field.signif_based_on_CV=1;
else
    signif_field.signif_based_on_CV=0;
end

%2. check if signif based on extreme bins:
benf_correction=alpha_val/length(r(1,:)); %correct for number of bins
benf_correction=benf_correction/2; %correct for checking both min and max

prct_alpha_max=prctile(r(2:end,:),100-benf_correction);
prct_alpha_min=prctile(r(2:end,:),benf_correction);
mid_prctile_for_width=prctile(r(2:end,:),width_at_heigth);
sig_bins=find(r(1,:)>prct_alpha_max | r(1,:)<prct_alpha_min);
if ~isempty(sig_bins)
    signif_field.signif_based_on_extreme_bins=1;
    signif_field.signif_bins=sig_bins;
    signif_field.pos_signif=find(r(1,:)>prct_alpha_max);
    signif_field.neg_signif=find(r(1,:)<prct_alpha_min);
    %find width of signif field at x height:
    %---------------------------------------
    length_data=length(r(1,:));
    data=r(1,:);
    [signif_field.width_pos,signif_field.width_line_x_pos,signif_field.width_line_y_pos,signif_field.pos_rise_time]=find_width(find(r(1,:)>mid_prctile_for_width),signif_field.pos_signif,data,dis_X_bins_vector_of_centers,mid_prctile_for_width);
    [signif_field.width_neg,signif_field.width_line_x_neg,signif_field.width_line_y_neg,signif_field.neg_rise_time]=find_width(find(r(1,:)<mid_prctile_for_width),signif_field.neg_signif,data,dis_X_bins_vector_of_centers,mid_prctile_for_width);
    % interpolate to get finer resolution:
    
    x_vec_interp=linspace(min(dis_X_bins_vector_of_centers),max(dis_X_bins_vector_of_centers),length(dis_X_bins_vector_of_centers)*10);
    length_data=length(x_vec_interp);
    non_nan_ind=~isnan(r(1,:));
    interpolate_r=interp1(dis_X_bins_vector_of_centers(non_nan_ind),r(1,non_nan_ind),x_vec_interp);
    data=interpolate_r;
    non_nan_ind=~isnan(mid_prctile_for_width);
    interpolate_mid=interp1(dis_X_bins_vector_of_centers(non_nan_ind),mid_prctile_for_width(non_nan_ind),x_vec_interp);
    pos_bins_val=dis_X_bins_vector_of_centers(signif_field.pos_signif);
    [~, interp_ind_pos]=min(abs(pos_bins_val-x_vec_interp'));
    [signif_field.width_pos_interp,signif_field.width_line_x_pos_interp,signif_field.width_line_y_pos_interp,signif_field.pos_rise_time_interp]=find_width(find(interpolate_r>interpolate_mid),interp_ind_pos,data,x_vec_interp,interpolate_mid);
    neg_bins_val=dis_X_bins_vector_of_centers(signif_field.neg_signif);
     [~, interp_ind_neg]=min(abs(neg_bins_val-x_vec_interp'));
    [signif_field.width_neg_interp,signif_field.width_line_x_neg_interp,signif_field.width_line_y_neg_interp,signif_field.neg_rise_time_interp]=find_width(find(interpolate_r<interpolate_mid),interp_ind_neg,data,x_vec_interp,interpolate_mid);
else
    signif_field.signif_based_on_extreme_bins=0;
end

signif_field.number_of_spikes_per_field=number_of_spikes;

%3. compute z scores:
zscore_shuffle=zscore(r);
signif_field.zscore_data=zscore_shuffle(1,:);

% save more data:

signif_field.shuffled_data=r(2:end,:);
signif_field.cv=cv(1);

    function [width,width_line_x_pos,width_line_y_pos,rise_time]=find_width(line_of_reference,ind_signif,data,x_bins,y_vec)
        length_data=length(data);
        non_nan_ind=find(~isnan(data));
        bin_size=mean(diff(x_bins));
        [ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(line_of_reference,length_data);
        if ~isempty(ind_length)
            for area_i=1:length(ind_length)
                %find if there is signif bin iind_signifn the area
                signif_bin_ind=find(ind_signif>start_ind(area_i) & ind_signif<end_ind(area_i));
                %if there is a signif field and it is not in the edge
                if ~isempty(signif_bin_ind) & start_ind(area_i)~=non_nan_ind(1) & end_ind(area_i)~=non_nan_ind(end)
                    width(area_i)=ind_length(area_i)*bin_size;
                    width_line_x_pos(area_i,:)=[x_bins(start_ind(area_i)) x_bins(end_ind(area_i))];
                    width_line_y_pos(area_i,:)=[y_vec(start_ind(area_i)) y_vec(end_ind(area_i))];
                    rise_time(area_i)=[x_bins(ind_signif(signif_bin_ind(1)))-x_bins(start_ind(area_i))];
                else
                    width(area_i)=nan;
                    width_line_x_pos(area_i,:)=[nan,nan];
                    width_line_y_pos(area_i,:)=[nan,nan];
                    rise_time(area_i)=nan;
                end
            end
        else
            width(1)=nan;
            width_line_x_pos(1,:)=[nan,nan];
            width_line_y_pos(1,:)=[nan,nan];
            rise_time(1)=nan;
        end
    end
end