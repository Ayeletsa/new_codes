%% shuffling for per field analysis
%-------------------------------------------------------------------------
%here we are doing non rigid shuffling, meaning that the spikes are not in
%the same temporal relation in each shuffle.
% we just take all spikes and randomly distribute them between all
% positions

function [signif_field]=shuffling_per_field_analysis_new_smooth(bsp,shuffles,dis_X_bins_vector,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D, frames_per_second,num_shuffles,alpha_val,old_smooth,smooth_window,smooth_type,smooth_tol,width_at_heigth,min_dis_pos_neg,distance,benf_correct);


% %behav_length=length(bsp_vec);
% parfor shuffle_i=1:num_shuffles
%     if shuffle_i==1
%         spike_vec_shuf=spike_mat;spike_vec_shuf=spike_vec_shuf(~isnan(spike_vec_shuf)); %original data
%         bsp_vec=bsp_mat;bsp_vec=bsp_vec(~isnan(bsp_vec)); %original data
%     else
%         perm_flights=randperm(size(bsp_mat,1));
%         [~, ind]=max(~isnan(bsp_mat),[],2);
%         intial_pos=bsp_mat(sub2ind(size(bsp_mat),1:size(bsp_mat,1),ind'));
%         initial_pos_permute=intial_pos(perm_flights);
%         intial_pos_to_move=initial_pos_permute'-intial_pos';
%         bsp_co_new=bsxfun(@plus,bsp_mat,intial_pos_to_move);
%         spike_co_new=bsxfun(@plus,spike_mat,intial_pos_to_move);
%         spike_vec_shuf=spike_co_new(:);spike_vec_shuf=spike_vec_shuf(~isnan(spike_vec_shuf));
%         bsp_vec=bsp_co_new(:);  bsp_vec=bsp_vec(~isnan(bsp_vec));
%   
% %         parfor co_i=1:length(perm_flights)
% %             
% %             bsp_dis_x_co=bsp_dis_mat(co_i,~isnan(bsp_dis_mat(co_i,:)));
% %             spike_dis_x_co=spike_dis_mat(co_i,~isnan(spike_dis_mat(co_i,:)));
% %             
% %             move_flight_to=bsp_dis_mat(perm_flights(co_i),~isnan(bsp_dis_mat(perm_flights(co_i),:)));
% %             %move_ego_by=median(move_flight_to)-median(bsp_dis_x_co);
% %             move_ego_by=(move_flight_to(1))-(bsp_dis_x_co(1));
% %             
% %             new_bsp_co{co_i}=bsp_dis_x_co+move_ego_by;
% %             new_spike_co{co_i}=spike_dis_x_co+move_ego_by;
% %             
% %             
% %         end
%         %bsp_vec=[new_bsp_co{:}];
%         %spike_vec_shuf=[new_spike_co{:}];
%         %         %choosing random ind for the spikes
%         %         rand_pos=randi(behav_length,number_of_spikes,1);
%         %         %distribute the spikes new pos
%         %         spike_vec_shuf=bsp_vec(rand_pos);
%         
%     end
if distance==1
    bsp_vec=bsp.dis_m;bsp_vec=bsp_vec(~isnan(bsp_vec));
else
    bsp_vec=bsp.time_to_co;bsp_vec=bsp_vec(~isnan(bsp_vec));
end
if size(bsp_vec,1)<size(bsp_vec,2)
    bsp_vec=bsp_vec';
end
%parfor shuffle_i=1:num_shuffles+1
    if distance==1
        spike_vec_shuf=[shuffles.dis_m];
    else
        spike_vec_shuf=[shuffles.time_to_co];
    end
    %compute tuning curve
     [~, ~, ~, r] ...
            = fn_compute_generic_1D_tuning_mat ...
            (bsp_vec,spike_vec_shuf',dis_X_bins_vector, time_spent_minimum_for_1D, frames_per_second);
       
%     [~, ~, ~, r(shuffle_i,:), ~,~] ...
%         = fn_compute_generic_1D_tuning_new_smooth ...
%         (bsp_vec,spike_vec_shuf,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D, frames_per_second, 0,0,0);
%     
    std_r=nanstd(r,[],2);
    mean_r=nanmean(r,2);
    cv=std_r./mean_r;
  
%end
spike_vec_shuf=shuffles(1).dis_m;spike_vec_shuf=spike_vec_shuf(~isnan(spike_vec_shuf));

number_of_spikes=sum(sum(~isnan(spike_vec_shuf)));

%% check signif
%1.check if signif based on cv:
prct_alpha=prctile(cv,100-alpha_val);
if cv(1)>prct_alpha
    signif_field.signif_based_on_CV=1;
else
    signif_field.signif_based_on_CV=0;
end

%2. check if signif based on extreme bins:
if benf_correct~=0
   benf_correction= alpha_val/benf_correct;
else
benf_correction=alpha_val/length(r(1,:)); %correct for number of bins
benf_correction=benf_correction/2; %correct for checking both min and max
end
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
    [signif_field.width_pos_interp,signif_field.width_line_x_pos_interp,signif_field.width_line_y_pos_interp,signif_field.pos_rise_time_interp,line_of_tunind_width_ind_pos]=find_width(find(interpolate_r>interpolate_mid),interp_ind_pos,data,x_vec_interp,interpolate_mid);
    neg_bins_val=dis_X_bins_vector_of_centers(signif_field.neg_signif);
    [~, interp_ind_neg]=min(abs(neg_bins_val-x_vec_interp'));
    [signif_field.width_neg_interp,signif_field.width_line_x_neg_interp,signif_field.width_line_y_neg_interp,signif_field.neg_rise_time_interp,line_of_tunind_width_ind_neg]=find_width(find(interpolate_r<interpolate_mid),interp_ind_neg,data,x_vec_interp,interpolate_mid);
else
    signif_field.signif_based_on_extreme_bins=0;
end

signif_field.number_of_spikes_per_field=number_of_spikes;

%3. compute z scores:
zscore_shuffle=zscore(r);
signif_field.zscore_data=zscore_shuffle(1,:);
%% check if coumpound (both neg and pos)
relevant_switch_times=[];
width_compund_all=[];
if ~distance
    if ~isempty(sig_bins)
        %find if there are close width
        %1. sum both pos and neg:
        bin_size_interp=min(diff(x_vec_interp));
        width_diff_time=diff(find(line_of_tunind_width_ind_pos+line_of_tunind_width_ind_neg)).*bin_size_interp;
        if any(width_diff_time<min_dis_pos_neg & width_diff_time>bin_size_interp)
            signif_field.compound=1;
            
            
        end
        if abs(min(signif_field.width_line_x_pos_interp(:))-max(signif_field.width_line_x_neg_interp(:)))<=min_dis_pos_neg | abs(max(signif_field.width_line_x_pos_interp(:))-min(signif_field.width_line_x_neg_interp(:)))<=min_dis_pos_neg
            signif_field.compound=1;
            %      compute switch time (time between last sig bin in pos to first in neg (or the other way)):
            if distance==0 % run only if time
                for pos_area_i=1:size(signif_field.width_line_x_pos_interp,1)
                    pos_position=signif_field.width_line_x_pos_interp(pos_area_i,:);
                    %find neg pos
                    neg_positions_all=signif_field.width_line_x_neg_interp;
                    for neg_pos_i=1:size(neg_positions_all,1)
                        if abs(min(pos_position)-max(neg_positions_all(neg_pos_i,:)))<=min_dis_pos_neg | abs(max(pos_position)-min(neg_positions_all(neg_pos_i,:)))<=min_dis_pos_neg
                            pos_ind_sig_bins_in_area=find(pos_bins_val>=pos_position(1) & pos_bins_val<=pos_position(2));
                            [~, interp_ind_pos]=min(abs(pos_bins_val(pos_ind_sig_bins_in_area)-x_vec_interp'));
                            
                            neg_ind_sig_bins_in_area=find(neg_bins_val>=neg_positions_all(neg_pos_i,1) & neg_bins_val<=neg_positions_all(neg_pos_i,2));
                            [~, interp_ind_neg]=min(abs(neg_bins_val(neg_ind_sig_bins_in_area)-x_vec_interp'));
                            
                            switch_times=[min(interp_ind_pos)-max(interp_ind_neg), min(interp_ind_neg)-max(interp_ind_pos)];
                            bin_size_for_switch=min(diff(x_vec_interp));
                            relevant_switch_times=[relevant_switch_times, (max(switch_times)*bin_size_for_switch)];
                            width_compund=signif_field.width_neg_interp(neg_pos_i)+signif_field.width_pos_interp(pos_area_i);
                            width_compund_all=[width_compund_all,width_compund];
                        end
                    end
                end
                signif_field.relevant_switch_times=relevant_switch_times;
                signif_field.compund_width=sum(width_compund_all);
            end
        else
            signif_field.compound=0;
        end
    else
        signif_field.compound=0;
    end
end
%% save more data:

signif_field.shuffled_data=r(2:end,:);
signif_field.cv=cv(1);
%%
    function [width,width_line_x_pos,width_line_y_pos,rise_time,line_of_tunind_width_ind]=find_width(line_of_reference,ind_signif,data,x_bins,y_vec)
        length_data=length(data);
        non_nan_ind=find(~isnan(data));
        bin_size=mean(diff(x_bins));
        [ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(line_of_reference,length_data);
        line_of_tunind_width_ind=zeros(1,length_data);
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
                    line_of_tunind_width_ind(start_ind(area_i):end_ind(area_i))=area_i;
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