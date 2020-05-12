function co = clicks_in_co(bsp_proc_data,behavioral_modes,tag_i,co_param_file_name,clicks_struct,plot_co_flag,audio_filt_struct,fig_folder,fig_prefix,audio_param_file_name,us2fs)
load(co_param_file_name)
load(audio_param_file_name)
us_factor=1e6;
%% create vectors of variables
self_ind = strcmp({clicks_struct.bat},'self');

% both bats are flying:
FE_self=bsp_proc_data(tag_i).flight_ind;
FE_other=bsp_proc_data(3-tag_i).flight_ind;
ind_both_bats_flying=intersect(FE_self,FE_other);
% find ind that bats are flying in different directions
temp_bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying,1)';
temp_bsp_other_x_pos = bsp_proc_data(3-tag_i).pos(ind_both_bats_flying,1)';
dir_flight=sign(diff(temp_bsp_x_pos));
dir_flight_other=sign(diff(temp_bsp_other_x_pos));
ind_different_direction=find((dir_flight-dir_flight_other)~=0)+1;
ind_both_bats_flying_different_dirs=ind_both_bats_flying(ind_different_direction);
    
% a. time stamps
bsp_ts = bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs);
clicks_ts = [clicks_struct(self_ind).clusters.peak_ts];

if plot_co_flag
    audio_ts = audio_filt_struct(1).ts;
    f=figure('WindowState','maximized');
end

% b. x position
bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1)';
clicks_x_pos = interp1(bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs) , bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1),clicks_ts);

% click intensity
clicks_intensity = [clicks_struct(self_ind).clusters.snr];

% double clicks
clicks_not_DC_flag = [clicks_struct(self_ind).clusters.not_DC];
ici = diff([clicks_struct(self_ind).clusters.peak_abs_idx]);
dc_1st_ind = find(min_intra_click_diff*us2fs < ici & ici < max_intra_click_diff*us2fs);
DC_counted_once = setdiff(1:length(clicks_ts),dc_1st_ind+1);

% relative distance from the other bat
bsp_dis_m = bsp_x_pos - bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1)';

% cross overs parameters
co_x_positions = bsp_proc_data(tag_i).pos(behavioral_modes.CO_point,1);
co_times_usec = bsp_proc_data(tag_i).ts(behavioral_modes.CO_point);
co_directions = sign(co_x_positions - bsp_proc_data(tag_i).pos(behavioral_modes.CO_point-1,1));
direction_ind = {find(co_directions>0),find(co_directions<0)};

% extract special clicks to plot
click_to_plot = struct('not_DC',[],'mode_2',[]);
for ibat = 1:2
    click_to_plot(ibat).not_DC = clicks_struct(ibat).clusters([clicks_struct(ibat).clusters.not_DC]>0);
    click_to_plot(ibat).mode_2 = clicks_struct(ibat).clusters([clicks_struct(ibat).clusters.mode2]>0);
end

%% loop for the 2 directions
for ii_dir = 1:2
    dir_ind = direction_ind{ii_dir};
    n_co_points = length(dir_ind);
    
    % create empty arrays for CO loop
    % a. bsp
    bsp_ts_usec = cell(n_co_points,1);
    bsp_time_to_co = cell(n_co_points,1);
    bsp_dis_m_at_co = cell(n_co_points,1);
    bsp_x_pos_at_co = cell(n_co_points,1);
    
    % b. clicks
    clicks_ts_usec = cell(n_co_points,1);
    clicks_time_to_co = cell(n_co_points,1);
    clicks_dis_m_at_co = cell(n_co_points,1);
    clicks_x_pos_at_co = cell(n_co_points,1);
    clicks_ind_from_total = cell(n_co_points,1);
    clicks_intensity_at_co = cell(n_co_points,1);
    clicks_not_DC_flag_co = cell(n_co_points,1);
    
    clicks_ts_usec_once = cell(n_co_points,1);
    clicks_time_to_co_once = cell(n_co_points,1);
    clicks_dis_m_at_co_once = cell(n_co_points,1);
    clicks_x_pos_at_co_once = cell(n_co_points,1);
    
    % loop for every cross over
    for ii_co = 1:n_co_points
        
        % find co time
        co_time_usec = co_times_usec(dir_ind(ii_co));
        
        % find bsp samples within a window of -+time_before_after from CO
        % a. find samples index
        bsp_relative_times = bsp_ts - co_time_usec;
        bsp_time_criteria_ind = ( - time_before_after_co*us_factor < bsp_relative_times & bsp_relative_times < time_before_after_co*us_factor);
        bsp_dis_criteria_ind = ( abs(bsp_dis_m)<dis_before_after_co)';
        bsp_ind = bsp_time_criteria_ind & bsp_dis_criteria_ind;
        % b. find samples values
        bsp_ts_usec{ii_co} =  bsp_ts(bsp_ind)';
        bsp_time_to_co{ii_co} = bsp_relative_times(bsp_ind)';
        bsp_dis_m_at_co{ii_co} = bsp_dis_m(bsp_ind);
        bsp_x_pos_at_co{ii_co} = bsp_x_pos(bsp_ind);
        
        % find clicks within a window of -+time_before_after from CO
        % a. find clicks index
        clicks_relative_times = clicks_ts - co_time_usec;
        clicks_ind = find(clicks_ts> bsp_ts_usec{ii_co}(1) & clicks_ts<bsp_ts_usec{ii_co}(end));

        % b. find clicks values
        clicks_ts_usec{ii_co} = clicks_ts(clicks_ind);
        clicks_time_to_co{ii_co} = clicks_relative_times(clicks_ind);
        clicks_dis_m_at_co{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_dis_m_at_co{ii_co},clicks_ts_usec{ii_co});
        clicks_x_pos_at_co{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_x_pos_at_co{ii_co},clicks_ts_usec{ii_co});
        clicks_intensity_at_co{ii_co} = clicks_intensity(clicks_ind);
        clicks_not_DC_flag_co{ii_co} = clicks_not_DC_flag(clicks_ind);
        
        % count double-clicks only once
        clicks_ind_once = intersect(clicks_ind,DC_counted_once);
        clicks_ts_usec_once{ii_co} = clicks_ts(clicks_ind_once);
        clicks_time_to_co_once{ii_co} = clicks_relative_times(clicks_ind_once);
        clicks_dis_m_at_co_once{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_dis_m_at_co{ii_co},clicks_ts_usec_once{ii_co});
        clicks_x_pos_at_co_once{ii_co} = interp1(bsp_ts_usec{ii_co},bsp_x_pos_at_co{ii_co},clicks_ts_usec_once{ii_co});
        
        % if direction is negative, flip relative positions so negative locations will represent positions before the CO
        if ii_dir == 2
            clicks_dis_m_at_co{ii_co} = - clicks_dis_m_at_co{ii_co};
            bsp_dis_m_at_co{ii_co} = - bsp_dis_m_at_co{ii_co};
            clicks_dis_m_at_co_once{ii_co} = - clicks_dis_m_at_co_once{ii_co};
        end
        
        if plot_co_flag
            C = {'g','r'};
            L = {'self bat logger (implanted)','other bat logger (not implanted)'};
            labels = {'ok','50 < snr < 70','70 < snr < 100','wrong length','wrong rise time','wrong spectrum','min diff','not DC','mode2 (lower spectrum)','close: snr < 100','close: aligned with other bat click'};
            
            %% plot audio signal for this CO
            co_ind_audio = audio_ts >= bsp_ts_usec{ii_co}(1) & audio_ts <= bsp_ts_usec{ii_co}(end);
            clf
            for ibat=1:2
                subplot(2,1,ibat)
                hold on
                co_ind_close = ismember([clicks_struct(ibat).clusters.peak_ts],audio_ts(co_ind_audio));
                co_ind_50 = co_ind_close & [clicks_struct(ibat).clusters.snr] < 70; 
                co_ind_70 = co_ind_close & [clicks_struct(ibat).clusters.snr] < 100 & [clicks_struct(ibat).clusters.snr] > 70;
                wrong_size = ismember([clicks_struct(ibat).wrong_size.peak_ts],audio_ts(co_ind_audio));
                wrong_rise = ismember([clicks_struct(ibat).wrong_rise_time.peak_ts],audio_ts(co_ind_audio));
                wrong_spect = ismember([clicks_struct(ibat).wrong_mid_low_spectrum.peak_ts],audio_ts(co_ind_audio));
                min_diff = ismember([clicks_struct(ibat).min_diff_reject.peak_ts],audio_ts(co_ind_audio));
                wrong_other_low = ismember([clicks_struct(ibat).wrong_other_bat_low.peak_ts],audio_ts(co_ind_audio));
                wrong_other_aligned = ismember([clicks_struct(ibat).wrong_other_bat_aligned.peak_ts],audio_ts(co_ind_audio));
                
                not_DC = ismember([click_to_plot(ibat).not_DC.peak_ts],audio_ts(co_ind_audio));
                mode_2 = ismember([click_to_plot(ibat).mode_2.peak_ts],audio_ts(co_ind_audio));
                
                plot((audio_ts(co_ind_audio) - co_time_usec)*1e-6,audio_filt_struct(ibat).filt(co_ind_audio)./clicks_struct(ibat).aud_BL,C{ibat});
                plot_co{1}=plot(([clicks_struct(ibat).clusters(co_ind_close).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).clusters(co_ind_close).snr],300),'k*');
                plot_co{2}=plot(([clicks_struct(ibat).clusters(co_ind_50).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).clusters(co_ind_50).snr],300),'m*');
                plot_co{3}=plot(([clicks_struct(ibat).clusters(co_ind_70).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).clusters(co_ind_70).snr],300),'c*');
                plot_co{4}=plot(([clicks_struct(ibat).wrong_size(wrong_size).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_size(wrong_size).snr],300),'bx');
                plot_co{5}=plot(([clicks_struct(ibat).wrong_rise_time(wrong_rise).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_rise_time(wrong_rise).snr],300),'b+');
                plot_co{6}=plot(([clicks_struct(ibat).wrong_mid_low_spectrum(wrong_spect).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_mid_low_spectrum(wrong_spect).snr],300),'b^');
                plot_co{7}=plot(([clicks_struct(ibat).min_diff_reject(min_diff).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).min_diff_reject(min_diff).snr],300),'bo');
                
                plot_co{8}=plot(([click_to_plot(ibat).not_DC(not_DC).peak_ts] - co_time_usec)*1e-6,min([click_to_plot(ibat).not_DC(not_DC).snr],300),'ko','MarkerSize',10);
                plot_co{9}=plot(([click_to_plot(ibat).mode_2(mode_2).peak_ts] - co_time_usec)*1e-6,min([click_to_plot(ibat).mode_2(mode_2).snr],300),'ks','MarkerSize',10);
                
                plot_co{10}=plot(([clicks_struct(ibat).wrong_other_bat_low(wrong_other_low).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_other_bat_low(wrong_other_low).snr],300),'yp');
                plot_co{11}=plot(([clicks_struct(ibat).wrong_other_bat_aligned(wrong_other_aligned).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_other_bat_aligned(wrong_other_aligned).snr],300),'y<');
                
                yline(th(ibat),'--');
                yline(-th(ibat),'--');
                ylim([-300,300])
                ylabel('SNR')
                xlabel('time relative to CO (s)')
                title(L{ibat})
                do_not_plot = cellfun('isempty',plot_co);
                legend([plot_co{~do_not_plot}],labels(~do_not_plot),'location','eastoutside')
            end
            fig_name = [fig_prefix '_co_dir_' num2str(ii_dir) '_co_' num2str(ii_co)];
            sgtitle(fig_name,'interpreter','none')
            saveas(f,fullfile(fig_folder,fig_name),'png')
            xlim(f.Children([3,5]),[-1 0.3])
            saveas(f,fullfile(fig_folder,[fig_name '_zoom']),'png')
            
            
            %% plot self-other alignment for this CO
            dis_criterion_for_click_alignment = bsp_dis_m_at_co{ii_co} >= click_offset_window(1) & bsp_dis_m_at_co{ii_co} <= click_offset_window(2);
            bsp_ts_click_offset = bsp_ts_usec{ii_co}(dis_criterion_for_click_alignment);
            co_ind_audio = audio_ts >= bsp_ts_click_offset(1) & audio_ts <= bsp_ts_click_offset(end);
            aud_ts_to_plot = (audio_ts(co_ind_audio) - co_time_usec)*1e-6;
            search_neighborhood = round(click_offset_search_window*us2fs);
            mean_snr = cellfun(@mean,{[clicks_struct(1).clusters.snr],[clicks_struct(2).clusters.snr]});
            
            clf
            ax{1} = axes('position',[0.04 0.5 0.95 0.35]);
            ax{2} = axes('position',[0.04 0.14 0.95 0.35]);
            ax{1}.XAxis.Visible = 'off';
            hold(ax{1},'on')
            hold(ax{2},'on')
            for ibat=1:2
                co_ind_close = ismember([clicks_struct(ibat).clicks_in_close_proximity.peak_ts],audio_ts(co_ind_audio));
                co_ind_all_clicks = ismember([clicks_struct(ibat).clusters.peak_ts],audio_ts(co_ind_audio));
                close_clicks_co = clicks_struct(ibat).clicks_in_close_proximity(co_ind_close);
                all_clicks_co = clicks_struct(ibat).clusters(co_ind_all_clicks);
                
                plot_aligned=plot(aud_ts_to_plot,audio_filt_struct(ibat).filt(co_ind_audio)./clicks_struct(ibat).aud_BL,C{ibat},'Parent',ax{ibat});
                plot(([all_clicks_co.peak_ts] - co_time_usec)*1e-6,min([all_clicks_co.snr],500),'k*','Parent',ax{ibat});
                yline(th(ibat),'--k','Parent',ax{ibat});
                yline(-th(ibat),'--k','Parent',ax{ibat});
                
                for cl=1:length(close_clicks_co)
                    idx_search = close_clicks_co(cl).abs_idx_at_2nd_bat;
                    search_start = (audio_ts(idx_search-search_neighborhood) - co_time_usec)*1e-6;
                    search_end = (audio_ts(idx_search+search_neighborhood) - co_time_usec)*1e-6;
                    search_amplitude = close_clicks_co(cl).snr / clicks_struct(ibat).mean_snr * clicks_struct(3-ibat).mean_snr;
%                     search_amplitude = th_for_aligned_clicks;
                    rectangle('Position',[search_start,-search_amplitude,search_end-search_start,2*search_amplitude],...
                        'FaceColor',[0.5 0.5 0.5 0.2],'EdgeColor','none','Parent',ax{3-ibat});
                end
                
                co_ind_wrong_aligned = ismember([clicks_struct(ibat).wrong_other_bat_aligned.peak_ts],audio_ts(co_ind_audio));
                co_ind_wrong_low = ismember([clicks_struct(ibat).wrong_other_bat_low.peak_ts],audio_ts(co_ind_audio));
                plot(([clicks_struct(ibat).wrong_other_bat_aligned(co_ind_wrong_aligned).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_other_bat_aligned(co_ind_wrong_aligned).snr],500),'b*','Parent',ax{ibat});
                plot(([clicks_struct(ibat).wrong_other_bat_low(co_ind_wrong_low).peak_ts] - co_time_usec)*1e-6,min([clicks_struct(ibat).wrong_other_bat_low(co_ind_wrong_low).snr],500),'bo','Parent',ax{ibat});
                
                ylim(ax{ibat},[-500 500])
                xlim(ax{ibat},[aud_ts_to_plot(1) aud_ts_to_plot(end)])
                ylabel('SNR')
                xlabel('time relative to CO (s)')
                legend(plot_aligned,L{ibat})
                uistack(plot_aligned,'top')
            end
            fig_name = [fig_prefix '_co_dir_' num2str(ii_dir) '_co_' num2str(ii_co) '_self_other_alignment'];
            title(fig_name,'interpreter','none','Parent',ax{1})
            saveas(f,fullfile(fig_folder,fig_name),'png')
        end
    end
    
    
    
    %% arrange data in matrices
    % bsp
    maxSize = max(cellfun(@numel,bsp_time_to_co));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
        % a. time to co
        rmat = cellfun(fcn,bsp_time_to_co,'UniformOutput',false);
        bsp.time_to_co = vertcat(rmat{:});
        % b. relative distance
        rmat = cellfun(fcn,bsp_dis_m_at_co,'UniformOutput',false);
        bsp.dis_m = vertcat(rmat{:});
    
    % clicks
    maxSize = max(cellfun(@numel,clicks_time_to_co));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
        % a. time to co
        % clicks_time_from_co(cellfun(@isempty,clicks_time_from_co)) = {nan};
        rmat = cellfun(fcn,clicks_time_to_co,'UniformOutput',false);
        clicks.time_to_co = vertcat(rmat{:});
        % b. relative distance
        % clicks_dis_m_at_co(cellfun(@isempty,clicks_dis_m_at_co)) = {nan};
        rmat = cellfun(fcn,clicks_dis_m_at_co,'UniformOutput',false);
        clicks.dis_m = vertcat(rmat{:});
        % c. intensity
        rmat = cellfun(fcn,clicks_intensity_at_co,'UniformOutput',false);
        clicks.intensity = vertcat(rmat{:});
        % single clicks
        rmat = cellfun(fcn,clicks_not_DC_flag_co,'UniformOutput',false);
        clicks.not_DC = vertcat(rmat{:});
    
    %% calculate calling rate
%     % 1. mean:
%         % a. time to co
%         bsp_vec_time = bsp.time_to_co(isfinite(bsp.time_to_co));
%         clicks_vec_time = clicks.time_to_co(isfinite(clicks.time_to_co));
%     
%         [~, ~, ~, time_to_co_fr, ~,~] ...
%             = fn_compute_generic_1D_tuning_new_smooth ...
%             (bsp_vec_time,clicks_vec_time,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
%     
%         calling_rate.time_to_co = [time_to_co_fr;time_X_bins_vector_of_centers];
%     
%         % b. relative distance
%         bsp_vec_dis = bsp.dis_m(isfinite(bsp.dis_m));
%         clicks_vec_dis = clicks.dis_m(isfinite(clicks.dis_m));
% 
%         [~, ~, ~, dis_m_fr, ~,~] ...
%             = fn_compute_generic_1D_tuning_new_smooth ...
%             (bsp_vec_dis,clicks_vec_dis,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
% 
%         calling_rate.dis_m = [dis_m_fr;dis_X_bins_vector_of_centers];
    
    % 2. single CO:
        % a. time to co
        compute_tuning_curve = @(x,y) fn_compute_generic_1D_tuning_new_smooth ...
            (x,y,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins/2, frames_per_second, 0,0,0);
    
        [~, ~, ~, single_co_time_to_co_fr, ~,~] = cellfun(compute_tuning_curve,bsp_time_to_co,clicks_time_to_co,'UniformOutput',false);
        calling_rate.single_co_time_to_co = cell2mat(single_co_time_to_co_fr);
        
        [~, ~, ~, single_co_time_to_co_fr_once, ~,~] = cellfun(compute_tuning_curve,bsp_time_to_co,clicks_time_to_co_once,'UniformOutput',false);
        calling_rate.single_co_time_to_co_once = cell2mat(single_co_time_to_co_fr_once);
        
        % b. relative distance
        compute_tuning_curve = @(x,y) fn_compute_generic_1D_tuning_new_smooth ...
            (x,y,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins/2, frames_per_second, 0,0,0);

        [~, ~, ~, single_co_dis_m_fr, ~,~] = cellfun(compute_tuning_curve,bsp_dis_m_at_co,clicks_dis_m_at_co,'UniformOutput',false);
        calling_rate.single_co_dis_m = cell2mat(single_co_dis_m_fr);
        
        [~, ~, ~, single_co_dis_m_fr_once, ~,~] = cellfun(compute_tuning_curve,bsp_dis_m_at_co,clicks_dis_m_at_co_once,'UniformOutput',false);
        calling_rate.single_co_dis_m_once = cell2mat(single_co_dis_m_fr_once);
    
    %% calculate click intensity
    % 1. mean:
%         intensity_vec = clicks.intensity(isfinite(clicks.intensity));
%         % a. time to co
%         calling_rate.intensity_time_to_co = tuning_curve_by_value_not_rate ...
%             (clicks_vec_time,bsp_vec_time,time_X_bins_vector,intensity_vec,time_spent_minimum_for_1D_bins,frames_per_second);
    
        % b. relative distance
%         calling_rate.intensity_dis_m = tuning_curve_by_value_not_rate ...
%             (clicks_vec_dis,bsp_vec_dis,amp_dis_X_bins_vector,intensity_vec,time_spent_minimum_for_1D_bins,frames_per_second,amp_binning_method);
%         
%         calling_rate.intensity_dis_m_bins = amp_dis_X_bins_vector_of_centers;
        
    % 2. single CO:
%         % a. time to co
%         fn_tuning_curve_val = @(x,y,z) tuning_curve_by_value_not_rate ...
%             (x,y,time_X_bins_vector,z,time_spent_minimum_for_1D_bins/2,frames_per_second);
%     
%         single_co_intensity_time_to_co = cellfun(fn_tuning_curve_val,clicks_time_to_co,bsp_time_to_co,clicks_intensity_at_co,'UniformOutput',false);
%         calling_rate.single_co_intensity_time_to_co = cell2mat(single_co_intensity_time_to_co);
    
        % b. relative distance
        fn_tuning_curve_val = @(x,y,z) tuning_curve_by_value_not_rate ...
            (x,y,amp_dis_X_bins_vector,z,time_spent_minimum_for_1D_bins/2,frames_per_second,amp_binning_method);

        single_co_intensity_dis_m = cellfun(fn_tuning_curve_val,clicks_dis_m_at_co,bsp_dis_m_at_co,clicks_intensity_at_co,'UniformOutput',false);
        calling_rate.single_co_intensity_dis_m = cell2mat(single_co_intensity_dis_m);
    
    
    %% save struct
%     co(ii_dir).clicks_ind = clicks_ind_from_total;
    co(ii_dir).clicks = clicks;
    co(ii_dir).bsp = bsp;
    co(ii_dir).calling_rate = calling_rate;
    
end
close
end