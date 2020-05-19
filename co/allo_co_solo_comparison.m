function solo_co_comparison=allo_co_solo_comparison(x_pos_bin_limits,sig_bins_width,solo_bsp_x_pos,solo_spikes_x_pos,co_bsp_x_pos,co_spikes_x_pos,frames_per_second,min_flights_per_bin)


    if ~isempty(x_pos_bin_limits)
        
     % compute bins:
     %----------------------------------------------------------
    x_pos_mod = mod(diff(x_pos_bin_limits),sig_bins_width);
    sig_bins_limits = x_pos_bin_limits(1)-x_pos_mod/2:sig_bins_width:x_pos_bin_limits(2)+x_pos_mod;
    n_bins = length(sig_bins_limits) - 1;
 
    % create empty variables
    p_fr = zeros(1,n_bins);
    diff_mean_fr = zeros(1,n_bins);
    n_solo_flights = size(solo_bsp_x_pos,1);
    n_co_points=size(co_bsp_x_pos,1);
    
    % Go over bins:
    for ii_bin = 1:n_bins
        
        bin_start = sig_bins_limits(ii_bin);
        bin_end = sig_bins_limits(ii_bin+1);
        solo_fr_at_bin = zeros(1,n_solo_flights);      
        co_fr_at_bin = zeros(1,n_co_points);
        
        %for every solo flight, calculate firing rate at bin
        %----------------------------------------------------------
        for ii_solo_flight = 1:n_solo_flights
            solo_bsp_at_flight = solo_bsp_x_pos(ii_solo_flight,:);
            n_bsp_at_bin = sum(bin_start <= solo_bsp_at_flight & solo_bsp_at_flight < bin_end);
            solo_spikes_at_flight = solo_spikes_x_pos(ii_solo_flight,:);
            n_spikes_at_bin = sum(bin_start <= solo_spikes_at_flight & solo_spikes_at_flight < bin_end);
            solo_fr_at_bin(ii_solo_flight) = (n_spikes_at_bin/n_bsp_at_bin) * frames_per_second;
        end
        
        %for every cross-over, calculate firing rate at bin
        %----------------------------------------------------------
        for ii_co_flight = 1:n_co_points
            co_bsp_at_flight = co_bsp_x_pos(ii_co_flight,:);
            n_bsp_at_bin = sum(bin_start <= co_bsp_at_flight & co_bsp_at_flight < bin_end);
            co_spikes_at_flight = co_spikes_x_pos(ii_co_flight,:);
            n_spikes_at_bin = sum(bin_start <= co_spikes_at_flight & co_spikes_at_flight < bin_end);
            co_fr_at_bin(ii_co_flight) = (n_spikes_at_bin/n_bsp_at_bin) * frames_per_second;
        end
        
        solo_fr_at_bin = solo_fr_at_bin(isfinite(solo_fr_at_bin));
        co_fr_at_bin = co_fr_at_bin(isfinite(co_fr_at_bin));
        
        % calculate p only for bins with more than n flights
        %----------------------------------------------------------
        if length(co_fr_at_bin) >= min_flights_per_bin
            [~,p_fr(ii_bin)] = ttest2(solo_fr_at_bin,co_fr_at_bin);
            diff_mean_fr(ii_bin) = mean(solo_fr_at_bin) - mean(co_fr_at_bin);
        else
            p_fr(ii_bin) = nan;
            diff_mean_fr(ii_bin) = nan;
        end
    end
    
    % save data to struct:
    %----------------------------------------------------------
    solo_co_comparison.p_diff_firing_rate = p_fr;
    solo_co_comparison.diff_mean_firing_rate = diff_mean_fr;
    sig_bins_centers = sig_bins_limits(1:end-1) + sig_bins_width/2;
    solo_co_comparison.bin_centers = sig_bins_centers;
    solo_co_comparison.n_bins = n_bins;

    
    else
       solo_co_comparison.n_bins=0; 
    end