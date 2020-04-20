function clicks = create_clicks_day_struct(clicks_struct_name,audio_filt_name,audio_param_file_name,bsp_proc_data,tag_i,behavioral_modes,co_param_file_name,solo_param_file_name,fig_folder)
%use audio data to detect clicks in solo and CO behaviors.
%use filtered audio!!!

load(audio_param_file_name)

[~,fig_prefix,~] = fileparts(clicks_struct_name);
fig_prefix = strrep(fig_prefix,'clicks_','');

% calculate distance between bats (for removing clicks from the other bat)
dist_self_from_other_x = (bsp_proc_data(tag_i).pos(:,1) - bsp_proc_data(3-tag_i).pos(:,1));
dist_self_from_other_x(behavioral_modes.directional_ind{2}) = -dist_self_from_other_x(behavioral_modes.directional_ind{2});

%% load audio and initialize struct
% clicks = struct();
% for ibat=1:2
%     if ibat==1
%             clicks(ibat).bat = 'self';
%             clicks(ibat).tag = tag_i;
%     else
%             clicks(ibat).bat = 'other';
%             clicks(ibat).tag = 3-tag_i;
%     end
%     
%     disp(['loading audio: ' clicks(ibat).bat])
%     fname = [audio_filt_name '_' clicks(ibat).bat '.ncs'];
%     [clicks(ibat).filt, clicks(ibat).ts, ~] = Nlx_csc_read(fname,[]);
%     
%     if strcmp(clicks(ibat).bat,'other') %use a single set of timestamps
%         clicks(ibat).filt = interp1(clicks(ibat).ts,clicks(ibat).filt,clicks(1).ts);
%         clicks(ibat).ts = clicks(1).ts;
%     end
%     
% end
% 
% % get fs_aud from the data
% rec_length=clicks(1).ts(end)-clicks(1).ts(1);
% unit=10^-(round(floor(log10(rec_length / 3600)) / 3)*3);
% fs_aud = 1/(unit*median(diff(clicks(1).ts(1:1000))));
% us2fs = fs_aud*1e-6; %convert usecs to samples for clarity

%% detect clicks
% I decided to detect all clicks, not only during COs
if ~exist(clicks_struct_name,'file')

for ibat=1:2
    disp(['detecting clicks: ' clicks(ibat).bat])
    
    % get flight epochs in audio (interpolate smaller gaps)
    intrp_wind = 3e6; %3secs
    clicks(ibat).ind_FE = interp1gap(bsp_proc_data(clicks(ibat).tag).ts_FE,ones(1,length(bsp_proc_data(clicks(ibat).tag).ts_FE)),clicks(ibat).ts,intrp_wind)>0;
    
    % calculate thresholds for detecting clicks (only during flights)
    nanind = isnan(clicks(ibat).filt);
    clicks(ibat).aud_BL = mad(clicks(ibat).filt(~nanind & clicks(ibat).ind_FE));
    
    %use simple threshold to find clicks
    thresh = th * clicks(ibat).aud_BL;
    filt_abs = abs(clicks(ibat).filt);
    above_th = filt_abs > thresh;
    
    %detect clicks only during flight times?
    above_th(~clicks(ibat).ind_FE) = 0;
    
    %connect points to get clusters
    intra_click_nulls = bwpropfilt(~above_th,'Area',[1 max_intra_click_null * us2fs]);
    above_th(intra_click_nulls) = 1;
    clicks(ibat).clusters = regionprops(above_th,filt_abs,'Area','PixelIdxList','PixelValues','MaxIntensity');
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
    
    %filter by amplitude
    snr_reject = [clicks(ibat).clusters.snr] < th;
    clicks(ibat).wrong_amplitude = clicks(ibat).clusters(snr_reject);
    clicks(ibat).clusters = clicks(ibat).clusters(~snr_reject);
    
    %identify double clicks
    ici = diff([clicks(ibat).clusters.peak_ts]);
    dc_1st_ind = find(min_intra_click_diff < ici & ici < max_intra_click_diff);
    dc_ind = union(dc_1st_ind,dc_1st_ind+1);
    not_dc = true(1,length(ici));
    not_dc(dc_ind) = 0;
    not_dc_c = num2cell([0 not_dc]);
    [clicks(ibat).clusters.not_DC] = not_dc_c{:};
    
    %calculate distance from 2nd bat
    dist_at_peak = interp1(bsp_proc_data(tag_i).ts,dist_self_from_other_x,[clicks(ibat).clusters.peak_ts]);
    %     if strcmp(clicks(i).bat,'other')
    %         dist_at_peak = -dist_at_peak;
    %     end
    dist_at_peak(dist_at_peak < click_offset_window(1) | click_offset_window(2) < dist_at_peak) = NaN;
    sound_offset = round(abs(dist_at_peak) / speed_of_sound * us2fs);
    abs_idx_at_2nd_bat = num2cell([clicks(ibat).clusters.peak_abs_idx] + sound_offset);
    [clicks(ibat).clusters.abs_idx_at_2nd_bat] = abs_idx_at_2nd_bat{:};
    
end

%% remove clicks originating from the other bat
search_neighborhood = round(click_offset_search_window*us2fs/2);
for ibat=1:2
    disp(['removing clicks on other bat: ' clicks(ibat).bat])
    
    clicks(ibat).clicks_in_close_proximity = clicks(ibat).clusters(~isnan([clicks(ibat).clusters.abs_idx_at_2nd_bat]));
    search_window = arrayfun(@(x) x-search_neighborhood:x+search_neighborhood, ...
        [clicks(ibat).clicks_in_close_proximity.abs_idx_at_2nd_bat],'UniformOutput',false);
    [~,suspected_ind_2nd_bat] = cellfun(@(x) ismember(x,[clicks(3-ibat).clusters.peak_abs_idx]),search_window,'UniformOutput',false);
    suspected_ind_2nd_bat = cellfun(@(x) x(x>0),suspected_ind_2nd_bat,'UniformOutput',false);
    suspected_ind_1st_bat = ~cellfun(@isempty,suspected_ind_2nd_bat);
    suspected_ind_2nd_bat = suspected_ind_2nd_bat(suspected_ind_1st_bat);
    
    clicks_with_lower_snr = cellfun(@(x,y) [clicks(3-ibat).clusters(x).snr] < y / clicks(ibat).mean_snr * clicks(3-ibat).mean_snr,...
        suspected_ind_2nd_bat,{clicks(ibat).clicks_in_close_proximity(suspected_ind_1st_bat).snr},'UniformOutput',false);
    has_wrong_click_on_2nd_bat = cell2mat(cellfun(@times,clicks_with_lower_snr,num2cell(1:length(suspected_ind_2nd_bat)),'UniformOutput',false));
    has_wrong_click_on_2nd_bat = has_wrong_click_on_2nd_bat(has_wrong_click_on_2nd_bat > 0);
    origin_clicks = clicks(ibat).clicks_in_close_proximity(suspected_ind_1st_bat);
    clicks(ibat).origin_clicks = origin_clicks(has_wrong_click_on_2nd_bat);
    
    suspected_ind_2nd_bat_vec = cell2mat(suspected_ind_2nd_bat);
    wrong_other_bat_clicks = suspected_ind_2nd_bat_vec(cell2mat(clicks_with_lower_snr));
    clicks(3-ibat).wrong_other_bat_clicks = clicks(3-ibat).clusters(wrong_other_bat_clicks);
end
% for i=1:2
%     clicks_to_remove = ismember([clicks(i).clusters.peak_ts],[clicks(i).wrong_other_bat_clicks.peak_ts]);
%     clicks(i).clusters = clicks(i).clusters(~clicks_to_remove);
% end

% %% plot clicks detected on other bat
% plot_wind = round(0.1e6*us2fs);
% figure('WindowState','maximized');
% C = {'g','r'};
% L = {'self bat logger (implanted)','other bat logger (not implanted)'};
% for i=1:2
%     disp(['plotting problematic clicks: ' clicks(i).bat])
%     for cl=1:length(clicks(3-i).wrong_other_bat_clicks)
%         idx_self = clicks(i).origin_clicks(cl).peak_abs_idx;
%         idx_search = clicks(i).origin_clicks(cl).abs_idx_at_2nd_bat;
%         idx_other = clicks(3-i).wrong_other_bat_clicks(cl).peak_abs_idx;
%         clf
%         hold on
% %         plot(clicks(i).ts(idx_self-plot_wind:idx_self+plot_wind),clicks(i).filt(idx_self-plot_wind:idx_self+plot_wind),C{i})
% %         plot(clicks(3-i).ts(idx_other-plot_wind:idx_other+plot_wind),clicks(3-i).filt(idx_other-plot_wind:idx_other+plot_wind),C{3-i})
%         plot(clicks(i).ts(idx_self-plot_wind:idx_self+plot_wind),clicks(i).filt(idx_self-plot_wind:idx_self+plot_wind)./clicks(i).aud_BL,C{i})
%         plot(clicks(3-i).ts(idx_other-plot_wind:idx_other+plot_wind),clicks(3-i).filt(idx_other-plot_wind:idx_other+plot_wind)./clicks(3-i).aud_BL,C{3-i})
% %         plot(clicks(3-i).wrong_other_bat_clicks(cl).peak_ts,clicks(3-i).wrong_other_bat_clicks(cl).MaxIntensity,'k*')
% %         plot(clicks(i).origin_clicks(cl).peak_ts,clicks(i).origin_clicks(cl).MaxIntensity,'b*')
%         plot(clicks(3-i).wrong_other_bat_clicks(cl).peak_ts,clicks(3-i).wrong_other_bat_clicks(cl).snr,'k*')
%         plot(clicks(i).origin_clicks(cl).peak_ts,clicks(i).origin_clicks(cl).snr,'b*')
%         search_start = clicks(3-i).ts(idx_search-search_neighborhood);
%         search_end = clicks(3-i).ts(idx_search+search_neighborhood);
% %         search_amplitude = clicks(i).origin_clicks(cl).snr*clicks(3-i).aud_BL;
%         search_amplitude = clicks(i).origin_clicks(cl).snr;
%         plot([clicks(3-i).ts(idx_search),clicks(3-i).ts(idx_search)],[-search_amplitude search_amplitude],'--b')
%         rectangle('Position',[search_start,-search_amplitude,search_end-search_start,2*search_amplitude],...
%             'FaceColor',[0 0 1 0.2],'EdgeColor','none');
%         fig_name = [fig_prefix '_' clicks(i).bat,'_wrong_click_',num2str(cl)];
% %         yline(-clicks(3-i).aud_BL*th,'--k');
% %         yline(clicks(3-i).aud_BL*th,'--k');
%         yline(th,'--k');
%         yline(-th,'--k');
%         legend(L([i 3-i]))
%         title({fig_name,['click originated from ' clicks(i).bat ', detected also on ' clicks(3-i).bat]},'Interpreter','none')
%         saveas(gcf,fullfile(fig_folder,fig_name),'png');
%     end
% end


%% save clicks struct
disp('saving clicks struct')
clicks_struct = rmfield(clicks,{'ts','filt','ind_FE'});
save(clicks_struct_name,'clicks_struct');

else
load(clicks_struct_name,'clicks_struct');
end

%% plot SNR histograms
figure
for ibat=1:2
    clicks_struct(ibat).mean_snr = mean([clicks_struct(ibat).clusters.snr]);
    subplot(1,2,ibat)
    histogram([clicks_struct(ibat).clusters.snr])
    xline(th,'r--');
    xline(clicks_struct(ibat).mean_snr,'k--');
    yl=ylim;
    text(clicks_struct(ibat).mean_snr + 10, yl(2) - 35,num2str(clicks_struct(ibat).mean_snr),'FontSize',15)
    title(clicks_struct(ibat).bat)
end
fig_name = [fig_prefix '_click_snr'];
sgtitle(fig_name,'interpreter','none')
saveas(gcf,fullfile(fig_folder,fig_name),'png');
close all

%% find clicks in solo & co + plot them
% %we do it here because the audio signals are already loaded to the memory
% clicks(1).co = clicks_in_co (bsp_proc_data,behavioral_modes,tag_i,co_param_file_name,clicks_struct,1,clicks,fig_folder,fig_prefix,audio_param_file_name,us2fs);
% clicks(1).solo = clicks_in_solo (bsp_proc_data,behavioral_modes,tag_i,solo_param_file_name,clicks_struct,1,clicks,fig_folder,fig_prefix,th);
