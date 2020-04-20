%load shuffles
%load initial co
%load co params

n_permutations = 1001;
up_factor = 10;

benf_correction=alpha_val/length(dis_X_bins_vector_of_centers); %correct for number of bins
benf_correction=benf_correction/2; %correct for checking both min and max

bsp = cell_co_solo_initial_analysis.co(ii_dir).bsp.dis_m;

spikes = cell_co_solo_initial_analysis.co(ii_dir).spikes.dis_m;
n_spikes_per_flight = sum(~isnan(spikes),2);
nflights = length(n_spikes_per_flight);

shuffled_spikes_dis = shuffling_struct(ii_dir).shuffled_data.spikes_dis_m;
n_shuffles = size(shuffled_spikes_dis,1);
shuffled_spikes_dis_flight = cell(1,nflights);

end_indices = cumsum(n_spikes_per_flight);
start_indices = end_indices - n_spikes_per_flight + 1;

for ii_flight = 1:nflights
    shuffled_spikes_dis_flight{ii_flight} = shuffled_spikes_dis(:,start_indices(ii_flight):end_indices(ii_flight));
end

%remove bad COs (high bl, no clear increment)
% shuffled_spikes_dis_flight = shuffled_spikes_dis_flight(~bad);
% nflights = length(shuffled_spikes_dis_flight);
% bsp(bad,:) = [];
co_early = randsample(nflights,floor(nflights/2));

perm_indices = zeros(floor(nflights/2),n_permutations);
delta = nan(1,n_permutations);
delta2 = nan(1,n_permutations);

parfor ii_perm = 2:n_permutations
    if ii_perm==1
        perm_ind_early = co_early;
    else
        perm_ind_early = randsample(nflights,floor(nflights/2));
    end
    perm_indices(:,ii_perm) = perm_ind_early;
    perm_ind_late = setdiff(1:nflights,perm_ind_early);
    
    ego_tuning_curve_early = calculate_rise_time_from_COs ... 
        (bsp,shuffled_spikes_dis_flight,perm_ind_early,benf_correction,up_factor, ...
        dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);

    ego_tuning_curve_late = calculate_rise_time_from_COs ... 
        (bsp,shuffled_spikes_dis_flight,perm_ind_late,benf_correction,up_factor, ...
        dis_X_bins_vector_of_centers,time_spent_minimum_for_1D_bins,frames_per_second);
    
    if ~isempty(ego_tuning_curve_early.switch_ind) && ~isempty(ego_tuning_curve_late.switch_ind)
        delta(ii_perm) = ego_tuning_curve_late.switch_ind - ego_tuning_curve_early.switch_ind;
    end
    if ~isempty(ego_tuning_curve_early.first_sig_bin) && ~isempty(ego_tuning_curve_late.first_sig_bin)
        delta2(ii_perm) = ego_tuning_curve_late.first_sig_bin - ego_tuning_curve_early.first_sig_bin;
    end
end

figure
histogram(delta)
figure
histogram(delta2)

function ego_tuning_curve = calculate_rise_time_from_COs(bsp,shuffled_spikes_flight,co2use,benf_correction,up_factor,x_bins,time_spent_minimum_for_1D,frames_per_second)
shuffled_spikes_dis = cell2mat(shuffled_spikes_flight(co2use));
n_shuffles = size(shuffled_spikes_dis,1);

bsp2use = bsp(co2use,:);
bsp_vec = bsp2use(~isnan(bsp2use(:)));

ego_shuf = nan(n_shuffles,length(x_bins));

for shuffle_i = 1:n_shuffles
    spike_vec_shuf = shuffled_spikes_dis(shuffle_i,:);

    [~, ~, ~, ego_shuf(shuffle_i,:), ~,~] = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec,spike_vec_shuf,x_bins, time_spent_minimum_for_1D, frames_per_second, 0,0,0);
end
prct_alpha_max=prctile(ego_shuf(2:end,:),100-benf_correction);
prct_alpha_min=prctile(ego_shuf(2:end,:),benf_correction);
mid_prctile_for_width=prctile(ego_shuf(2:end,:),50);
first_sig_bin=find(ego_shuf(1,:)>prct_alpha_max | ego_shuf(1,:)<prct_alpha_min,1);

if ~isempty(first_sig_bin)
    x_bins_up=linspace(min(x_bins),max(x_bins),length(x_bins)*up_factor);
    non_nan_ind=~isnan(ego_shuf(1,:));
    ego_true_up=interp1(x_bins(non_nan_ind),ego_shuf(1,non_nan_ind),x_bins_up);
    non_nan_ind=~isnan(mid_prctile_for_width);
    mid_shuf_up=interp1(x_bins(non_nan_ind),mid_prctile_for_width(non_nan_ind),x_bins_up);
    
    if any(ego_shuf(1,first_sig_bin)>prct_alpha_max)
        bin_type = 'pos';
        lower_than_mid = find(ego_true_up < mid_shuf_up);
        switch_ind = find(lower_than_mid < first_sig_bin,1,'last');
    elseif any(ego_shuf(1,first_sig_bin)<prct_alpha_min)
        bin_type = 'neg';
        higher_than_mid = find(ego_true_up > mid_shuf_up);
        switch_ind = find(higher_than_mid < first_sig_bin,1,'last');
    end
    
    ego_tuning_curve.first_sig_bin = first_sig_bin;
    ego_tuning_curve.sign = bin_type;
    ego_tuning_curve.switch_ind = switch_ind;
    
else
    ego_tuning_curve.first_sig_bin = [];
    ego_tuning_curve.switch_ind = [];
end
end

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
        if ~isempty(signif_bin_ind) && start_ind(area_i)~=non_nan_ind(1) && end_ind(area_i)~=non_nan_ind(end)
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

