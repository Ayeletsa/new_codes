function create_clicks_day_struct(p,clicks_struct_name,audio_filt_name,audio_param_file_name,general_behavior_data_file_name,tag_i,behav_struct_name,co_param_file_name,solo_param_file_name,fig_folder)
%use audio data to detect clicks in solo and CO behaviors.
%use filtered audio!!!

if exist(clicks_struct_name,'file')
    skip_detection = true;
else
    skip_detection = false;
end

load(audio_param_file_name)
load(behav_struct_name)
load(general_behavior_data_file_name,'bsp_proc_data')
    
[~,fig_prefix,~] = fileparts(clicks_struct_name);
% fig_prefix = strrep(fig_prefix,'clicks_',['null' num2str(max_intra_click_null) '_']);
fig_prefix = strrep(fig_prefix,'clicks_','');

%for problematic days we have slightly different parameters:
if ~contains(fig_prefix,problematic_days)
    th(1)=th(2);
    min_band_energy_ratio_low_snr(1)=min_band_energy_ratio_low_snr(2);
    min_band_energy_ratio_high_snr(1)=min_band_energy_ratio_high_snr(2);
    th_for_aligned_clicks(1)=th_for_aligned_clicks(2);
    th_for_close_clicks(1)=th_for_close_clicks(2);
end

% calculate distance between bats (for removing clicks from the other bat)
dist_self_from_other_x = (bsp_proc_data(tag_i).pos(:,1) - bsp_proc_data(3-tag_i).pos(:,1));
dist_self_from_other_x(behavioral_modes.directional_ind{2}) = -dist_self_from_other_x(behavioral_modes.directional_ind{2});

%% load audio and initialize struct
clicks = struct('name',fig_prefix,'audio_params',audio_params);
for ibat=1:2
    if ibat==1
            clicks(ibat).bat = 'self';
            clicks(ibat).tag = tag_i;
            unfilt_fname = fullfile(p.path_day_dir,p.Audio_dir_self,'audio.ncs');
    else
            clicks(ibat).bat = 'other';
            clicks(ibat).tag = 3-tag_i;
            unfilt_fname = fullfile(p.path_day_dir,p.Audio_dir_other,'audio.ncs');
    end
    
    disp(['loading audio: ' clicks(ibat).bat])
    fname = [audio_filt_name '_' clicks(ibat).bat '.ncs'];
    [clicks(ibat).filt, clicks(ibat).ts, ~] = Nlx_csc_read(fname,[]);
    if ~skip_detection
    [clicks(ibat).unfilt, ~, ~] = Nlx_csc_read(unfilt_fname,[]);
    end
    if strcmp(clicks(ibat).bat,'other') %use a single set of timestamps
        clicks(ibat).filt = interp1(clicks(ibat).ts,clicks(ibat).filt,clicks(1).ts);
        if ~skip_detection
        clicks(ibat).unfilt = interp1(clicks(ibat).ts,clicks(ibat).unfilt,clicks(1).ts);
        end
        clicks(ibat).ts = clicks(1).ts;
    end
end

% get fs_aud from the data
rec_length=clicks(1).ts(end)-clicks(1).ts(1);
unit=10^-(round(floor(log10(rec_length / 3600)) / 3)*3);
fs_aud = 1/(unit*median(diff(clicks(1).ts(1:1000))));
us2fs = fs_aud*1e-6; %convert usecs to samples for clarity

%% detect clicks
% I decided to detect all clicks, not only during COs
if ~skip_detection
for ibat=1:2
    disp(['detecting clicks: ' clicks(ibat).bat])
    
    % get flight epochs in audio (interpolate smaller gaps)
    intrp_wind = 3e6; %3secs
    clicks(ibat).ind_FE = interp1gap(bsp_proc_data(clicks(ibat).tag).ts_FE,ones(1,length(bsp_proc_data(clicks(ibat).tag).ts_FE)),clicks(ibat).ts,intrp_wind)>0;
    
    % calculate thresholds for detecting clicks (only during flights)
    nanind = isnan(clicks(ibat).filt);
    clicks(ibat).aud_BL = mad(clicks(ibat).filt(~nanind & clicks(ibat).ind_FE));
    
    %use simple threshold to find clicks
    thresh = th(ibat) * clicks(ibat).aud_BL;
    filt_abs = abs(clicks(ibat).filt);
    above_th = filt_abs > thresh;
    
    %detect clicks only during flight times?
    above_th(~clicks(ibat).ind_FE) = 0;
    
    %connect points to get clusters
    intra_click_nulls = bwpropfilt(~above_th,'Area',[1 max_intra_click_null * us2fs]);
    above_th(intra_click_nulls) = 1;
    clicks(ibat).clusters = regionprops(above_th,filt_abs,'Area','PixelIdxList','PixelValues','MaxIntensity');
    
    clear nanind above_th filt_abs intra_click_nulls
    
    snr = num2cell([clicks(ibat).clusters.MaxIntensity]/clicks(ibat).aud_BL);
    [clicks(ibat).clusters.snr] = snr{:};
    clicks(ibat).mean_snr = mean([clicks(ibat).clusters.snr]);
    
    [~,peak_idx]=cellfun(@max,{clicks(ibat).clusters.PixelValues},'UniformOutput',false);
    cluster_start = cellfun(@(x) x(1),{clicks(ibat).clusters.PixelIdxList},'UniformOutput',false);
    cluster_end = cellfun(@(x) x(end),{clicks(ibat).clusters.PixelIdxList},'UniformOutput',false);
    peak_abs_idx = cellfun(@(x,y) x+y-1 ,peak_idx,cluster_start,'UniformOutput',false);
    peak_ts = cellfun(@(x) clicks(ibat).ts(x),peak_abs_idx,'UniformOutput',false);
    [clicks(ibat).clusters.start] = cluster_start{:};
    [clicks(ibat).clusters.end] = cluster_end{:};
    [clicks(ibat).clusters.peak_abs_idx] = peak_abs_idx{:};
    [clicks(ibat).clusters.peak_ts] = peak_ts{:};
    
    %filter by rise times
    rise_time_reject = cell2mat(peak_idx) > max_rise_time * us2fs;
    clicks(ibat).wrong_rise_time = clicks(ibat).clusters(rise_time_reject);
    clicks(ibat).clusters = clicks(ibat).clusters(~rise_time_reject);
    
    %filter clusters by size
    size_reject = [clicks(ibat).clusters.Area] < min_cluster_size * us2fs | max_cluster_size * us2fs < [clicks(ibat).clusters.Area];
    clicks(ibat).wrong_size = clicks(ibat).clusters(size_reject);
    clicks(ibat).clusters = clicks(ibat).clusters(~size_reject);
    
%     %filter by amplitude
%     snr_reject = [clicks(ibat).clusters.snr] < th;
%     clicks(ibat).wrong_amplitude = clicks(ibat).clusters(snr_reject);
%     clicks(ibat).clusters = clicks(ibat).clusters(~snr_reject);
    
    %Shir's spectral criterion:
    %filter by frequency domain information    
    fft_window = (time_window_for_FFT(1)*us2fs): ...
    (time_window_for_FFT(2)*us2fs + 1);
    L = length(fft_window);
    for ii_click = 1:length(clicks(ibat).clusters)
        ind_fft = clicks(ibat).clusters(ii_click).peak_abs_idx + fft_window;
        fft_click = fft(clicks(ibat).unfilt(ind_fft));
        fft_click = abs(fft_click/L);
        fft_click_single_side = fft_click(1:L/2+1);
        fft_click_single_side(2:end-1) = 2*fft_click_single_side(2:end-1);
        f_fft = (fs_aud*(0:(L/2))/L)*10^-3; %KHz;
        low_freq_ind = and(f_fft>low_freq_band(1), ...
            f_fft<low_freq_band(2));
        mid_freq_ind = and(f_fft>mid_freq_band(1), ...
            f_fft<mid_freq_band(2));
        high_freq_ind = and(f_fft>high_freq_band(1), ...
            f_fft<high_freq_band(2));
        clicks(ibat).clusters(ii_click).spectrum_ratio_mid_low = ...
            mag2db(mean(fft_click_single_side(mid_freq_ind))/mean(fft_click_single_side(low_freq_ind)));
%         clicks(ibat).clusters(ii_click).spectrum_ratio_high_mid = ...
%             mag2db(mean(fft_click_single_side(high_freq_ind))/mean(fft_click_single_side(mid_freq_ind)));
    end
    %we use different criteria based on clicks snr
    mid_low_spectrum_reject_low_snr = [clicks(ibat).clusters.spectrum_ratio_mid_low] < min_band_energy_ratio_low_snr(ibat);
    mid_low_spectrum_reject_high_snr = [clicks(ibat).clusters.spectrum_ratio_mid_low] < min_band_energy_ratio_high_snr(ibat);
      
    high_snr_clicks = [clicks(ibat).clusters.snr] > snr_th_for_energy_ratio;
    wrong_mid_low_spectrum = (mid_low_spectrum_reject_high_snr & high_snr_clicks) | (mid_low_spectrum_reject_low_snr & ~high_snr_clicks);
    
    %identify "mode2/weird spectrum" clicks (high snr + lower spectrum)
    weird_spect = num2cell(high_snr_clicks & mid_low_spectrum_reject_low_snr & ~mid_low_spectrum_reject_high_snr);
    [clicks(ibat).clusters.mode2] = weird_spect{:};
    
    clicks(ibat).wrong_mid_low_spectrum = clicks(ibat).clusters(wrong_mid_low_spectrum);
    clicks(ibat).clusters = clicks(ibat).clusters(~wrong_mid_low_spectrum);
    
    %remove clicks too close to other clicks (echos or noises - chooses the
    % click with higher snr)
    M = zeros(1,length(clicks(ibat).filt));
    M([clicks(ibat).clusters.peak_abs_idx]) = [clicks(ibat).clusters.snr];
    [~,new_click_IX] = findpeaks(M, 'MinPeakDistance', min_click_diff*us2fs);
    [~, min_diff_reject] = setdiff([clicks(ibat).clusters.peak_abs_idx], new_click_IX);
    clicks(ibat).min_diff_reject = clicks(ibat).clusters(min_diff_reject);
    clicks(ibat).clusters(min_diff_reject) = [];

    %calculate distance from 2nd bat
    dist_at_peak = interp1(bsp_proc_data(tag_i).ts,dist_self_from_other_x,[clicks(ibat).clusters.peak_ts]);
    dist_at_peak(dist_at_peak < click_offset_window(1) | click_offset_window(2) < dist_at_peak) = NaN;
    sound_offset = round(abs(dist_at_peak) / speed_of_sound * us2fs);
    abs_idx_at_2nd_bat = num2cell([clicks(ibat).clusters.peak_abs_idx] + sound_offset);
    [clicks(ibat).clusters.abs_idx_at_2nd_bat] = abs_idx_at_2nd_bat{:};
    
end

%% remove clicks originating from the other bat
search_neighborhood = round(click_offset_search_window*us2fs);
for ibat=1:2
    disp(['removing clicks on other bat: ' clicks(ibat).bat])
    
    close_proximity_ind = ~isnan([clicks(ibat).clusters.abs_idx_at_2nd_bat]);
    clicks(ibat).clicks_in_close_proximity = clicks(ibat).clusters(close_proximity_ind);
    clicks(ibat).clicks_in_close_proximity([clicks(ibat).clicks_in_close_proximity.snr] < th_for_aligned_clicks(ibat)) = [];
    search_window = arrayfun(@(x) x-search_neighborhood:x+search_neighborhood, ...
        [clicks(ibat).clicks_in_close_proximity.abs_idx_at_2nd_bat],'UniformOutput',false);
%     [~,suspected_ind_2nd_bat] = cellfun(@(x) ismember(x,[clicks(3-ibat).clusters.peak_abs_idx]),search_window,'UniformOutput',false);
%     suspected_ind_2nd_bat = cellfun(@(x) x(x>0),suspected_ind_2nd_bat,'UniformOutput',false);
%     suspected_ind_1st_bat = ~cellfun(@isempty,suspected_ind_2nd_bat);
%     suspected_ind_2nd_bat = suspected_ind_2nd_bat(suspected_ind_1st_bat);
    
%     clicks_with_lower_snr = cellfun(@(x,y) [clicks(3-ibat).clusters(x).snr] < y / clicks(ibat).mean_snr * clicks(3-ibat).mean_snr,...
%         suspected_ind_2nd_bat,{clicks(ibat).clicks_in_close_proximity(suspected_ind_1st_bat).snr},'UniformOutput',false); 
%     has_wrong_click_on_2nd_bat = cell2mat(cellfun(@times,clicks_with_lower_snr,num2cell(1:length(suspected_ind_2nd_bat)),'UniformOutput',false));
%     has_wrong_click_on_2nd_bat = has_wrong_click_on_2nd_bat(has_wrong_click_on_2nd_bat > 0);
%     origin_clicks = clicks(ibat).clicks_in_close_proximity(suspected_ind_1st_bat);
%     clicks(ibat).origin_clicks = origin_clicks(has_wrong_click_on_2nd_bat);
%     
%     suspected_ind_2nd_bat_vec = cell2mat(suspected_ind_2nd_bat);
%     wrong_other_bat_clicks = suspected_ind_2nd_bat_vec(cell2mat(clicks_with_lower_snr));
%     clicks(3-ibat).wrong_other_bat_clicks = clicks(3-ibat).clusters(wrong_other_bat_clicks);

    all_search_windows = cell2mat(search_window);
    aligned_clicks_ind = ismember([clicks(3-ibat).clusters.peak_abs_idx],all_search_windows);
    clicks_with_low_snr_aligned = [clicks(3-ibat).clusters.snr] < th_for_aligned_clicks(3-ibat);
    clicks_with_low_snr_close = [clicks(ibat).clusters.snr] < th_for_close_clicks(ibat);
    clicks(3-ibat).wrong_other_bat_aligned = clicks(3-ibat).clusters(aligned_clicks_ind & clicks_with_low_snr_aligned);
    clicks(ibat).wrong_other_bat_low = clicks(ibat).clusters(close_proximity_ind & clicks_with_low_snr_close);  
end
for ibat=1:2
    clicks_to_remove = ismember([clicks(ibat).clusters.peak_ts],cat(2,[clicks(ibat).wrong_other_bat_aligned.peak_ts],[clicks(ibat).wrong_other_bat_low.peak_ts]));
    clicks(ibat).clusters = clicks(ibat).clusters(~clicks_to_remove);
end


%% identify double clicks
%this is sensitive to clicks from the other bat, so we do it afterwards
for ibat=1:2
    ici = diff([clicks(ibat).clusters.peak_abs_idx]);
    dc_1st_ind = find(min_intra_click_diff*us2fs < ici & ici < max_intra_click_diff*us2fs);
    dc_ind = union(dc_1st_ind,dc_1st_ind+1);
    not_dc = true(1,length(ici)+1);
    not_dc(dc_ind) = 0;
    not_dc_cell = num2cell(not_dc);
    [clicks(ibat).clusters.not_DC] = not_dc_cell{:};
end

%% save clicks struct
clicks_struct = rmfield(clicks,{'ts','filt','unfilt','ind_FE'});
save(clicks_struct_name,'clicks_struct');
else
disp('skip detection')    
load(clicks_struct_name,'clicks_struct')
end

%% plot SNR histograms
% figure('WindowState','maximized')
% C={'g','r'};
% hist_end = 500;
% bin_size = 5;
% for ibat=1:2
%     clicks_struct(ibat).mean_snr = mean([clicks_struct(ibat).clusters.snr]);
%     subplot(2,1,ibat)
%     histogram([clicks_struct(ibat).clusters.snr],'BinEdges',0:bin_size:hist_end,'FaceColor',C{ibat})
%     xline(th,'--');
% %     xline(clicks_struct(ibat).mean_snr,'k--');
% %     yl=ylim;
% %     text(clicks_struct(ibat).mean_snr + 10, yl(2) - 35,num2str(clicks_struct(ibat).mean_snr),'FontSize',15)
% %     text(th+10,yl(2)-35,num2str(th),'FontSize',15)
%     title(clicks_struct(ibat).bat)
%     xlim([th-bin_size,hist_end])
% end
% fig_name = [fig_prefix '_click_snr'];
% sgtitle([fig_name '_th_' num2str(th)],'interpreter','none')
% saveas(gcf,fullfile(fig_folder,fig_name),'png');
% close all

%% plot click_statistics
clicks_stats
% plot_click_spectrum

%% find clicks in solo & co + plot them
% %we do it here because the audio signals are already loaded to the memory
clicks_struct(1).co = clicks_in_co (bsp_proc_data,behavioral_modes,tag_i,co_param_file_name,clicks_struct,0,clicks,fig_folder,fig_prefix,audio_param_file_name,us2fs);
clicks_struct(1).solo = clicks_in_solo (bsp_proc_data,behavioral_modes,tag_i,solo_param_file_name,clicks_struct,1,clicks,fig_folder,fig_prefix,audio_param_file_name,us2fs);

save(clicks_struct_name,'clicks_struct');