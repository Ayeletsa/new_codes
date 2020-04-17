%% shuffling for per field analysis
%-------------------------------------------------------------------------
%here we are doing non rigid shuffling, meaning that the spikes are not in
%the same temporal relation in each shuffle.
% we just take all spikes and randomly distribute them between all
% positions

function [signif_field]=shuffling_per_field_analysis(bsp_vec,spikes_vec,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D, frames_per_second,num_shuffles,alpha_val);
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
        = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spike_vec_shuf,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D, frames_per_second, 0,0,0);
    
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

prct_alpha_max=prctile(r,100-benf_correction);
prct_alpha_min=prctile(r,benf_correction);

sig_bins=find(r(1,:)>prct_alpha_max | r(1,:)<prct_alpha_min);
if ~isempty(sig_bins)
    signif_field.signif_based_on_extreme_bins=1;
    signif_field.signif_bins=sig_bins;
     signif_field.pos_signif=find(r(1,:)>prct_alpha_max);
     signif_field.neg_signif=find(r(1,:)<prct_alpha_min);
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
