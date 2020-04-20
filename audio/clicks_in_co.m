function co = clicks_in_co(bsp_proc_data,behavioral_modes,tag_i,co_param_file_name,clicks_struct,plot_co_flag,audio_filt_struct,fig_folder,fig_prefix,audio_param_file_name,us2fs)
load(co_param_file_name)
us_factor=1e6;
%% create vectors of variables

% both bats are flying:
use_flight_criterion = 1;
if use_flight_criterion
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
    
else
    ind_both_bats_flying_different_dirs = true(1,length(bsp_proc_data(tag_i).ts));
end

% a. time stamps
bsp_ts = bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs);
clicks_ts = [clicks_struct(strcmp({clicks_struct.bat},'self')).clusters.peak_ts];
if plot_co_flag
    audio_ts = audio_filt_struct(1).ts;
    load(audio_param_file_name)
end

% b. x position
bsp_x_pos = bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1)';
clicks_x_pos = interp1(bsp_proc_data(tag_i).ts(ind_both_bats_flying_different_dirs) , bsp_proc_data(tag_i).pos(ind_both_bats_flying_different_dirs,1),clicks_ts);

% click intensity
clicks_intensity = [clicks_struct(strcmp({clicks_struct.bat},'self')).clusters.snr];

% relative distance from the other bat
bsp_dis_m = bsp_x_pos - bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1)';
other_bat_pos_at_clicks = interp1(bsp_proc_data(3-tag_i).ts(ind_both_bats_flying_different_dirs) , bsp_proc_data(3-tag_i).pos(ind_both_bats_flying_different_dirs,1),clicks_ts);
clicks_dis_m = clicks_x_pos - other_bat_pos_at_clicks;

% cross overs parameters
co_x_positions = bsp_proc_data(tag_i).pos(behavioral_modes.CO_point,1);
co_times_usec = bsp_proc_data(tag_i).ts(behavioral_modes.CO_point);
co_directions = sign(co_x_positions - bsp_proc_data(tag_i).pos(behavioral_modes.CO_point-1,1));
direction_ind = {find(co_directions>0),find(co_directions<0)};

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
    
    f=figure('WindowState','maximized');
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
        clicks_time_criteria_ind = ( - time_before_after_co*us_factor < clicks_relative_times & clicks_relative_times < time_before_after_co*us_factor);
        clicks_dis_criteria_ind = ( abs( clicks_dis_m)<dis_before_after_co);
        clicks_ind = clicks_time_criteria_ind & clicks_dis_criteria_ind;
        % b. find clicks values
        clicks_ts_usec{ii_co} = clicks_ts(clicks_ind);
        clicks_time_to_co{ii_co} = clicks_relative_times(clicks_ind);
        clicks_dis_m_at_co{ii_co} = clicks_dis_m(clicks_ind);
        clicks_x_pos_at_co{ii_co} = clicks_x_pos(clicks_ind);
        clicks_ind_from_total{ii_co} = find(clicks_ind);
        clicks_intensity_at_co{ii_co} = clicks_intensity(clicks_ind);
        
        % if direction is negative, flip relative positions so negative locations will represent positions before the CO
        if ii_dir == 2
            clicks_dis_m_at_co{ii_co} = - clicks_dis_m_at_co{ii_co};
            bsp_dis_m_at_co{ii_co} = - bsp_dis_m_at_co{ii_co};
        end
        
        if plot_co_flag
            %% plot audio signal for this CO
            co_ind_audio = audio_ts >= bsp_ts_usec{ii_co}(1) & audio_ts <= bsp_ts_usec{ii_co}(end);
            clf
            C = {'g','r'};
            L = {'self bat logger (implanted)','other bat logger (not implanted)'};
            for i=1:2
                subplot(2,1,i)
                hold on
                co_ind_click = ismember([clicks_struct(i).clusters.peak_ts],audio_ts(co_ind_audio));
                plot((audio_ts(co_ind_audio) - co_time_usec)*1e-6,audio_filt_struct(i).filt(co_ind_audio)./clicks_struct(i).aud_BL,C{i})
                p=plot(([clicks_struct(i).clusters(co_ind_click).peak_ts] - co_time_usec)*1e-6,[clicks_struct(i).clusters(co_ind_click).snr],'k*');
                yline(th,'--');
                yline(-th,'--');
                set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                ylabel('SNR')
                xlabel('time relative to CO (s)')
                title(L{i})
            end
            fig_name = [fig_prefix '_co_dir_' num2str(ii_dir) '_co_' num2str(ii_co)];
            sgtitle(fig_name,'interpreter','none')
            saveas(f,fullfile(fig_folder,fig_name),'png')
            
            %% plot self-other alignment for this CO
            dis_criterion_for_click_alignment = bsp_dis_m_at_co{ii_co} >= click_offset_window(1) & bsp_dis_m_at_co{ii_co} <= click_offset_window(2);
            bsp_ts_click_offset = bsp_ts_usec{ii_co}(dis_criterion_for_click_alignment);
            co_ind_audio = audio_ts >= bsp_ts_click_offset(1) & audio_ts <= bsp_ts_click_offset(end);
            search_neighborhood = round(click_offset_search_window*us2fs/2);
            mean_snr = cellfun(@mean,{[clicks_struct(1).clusters.snr],[clicks_struct(2).clusters.snr]});
            
            clf
            ax{1} = axes('position',[0.1 0.5 0.8 0.35]);
            ax{2} = axes('position',[0.1 0.14 0.8 0.35]);
            ax{1}.XAxis.Visible = 'off';
            hold(ax{1},'on')
            hold(ax{2},'on')
            for i=1:2
                co_ind_click = ismember([clicks_struct(i).clusters.peak_ts],audio_ts(co_ind_audio));
                close_clicks_co = clicks_struct(i).clusters(co_ind_click);
                p=plot((audio_ts(co_ind_audio) - co_time_usec)*1e-6,audio_filt_struct(i).filt(co_ind_audio)./clicks_struct(i).aud_BL,C{i},'Parent',ax{i});
                plot(([close_clicks_co.peak_ts] - co_time_usec)*1e-6,[close_clicks_co.snr],'k*','Parent',ax{i});
                yline(th,'--k','Parent',ax{i});
                yline(-th,'--k','Parent',ax{i});
                
                for cl=1:length(close_clicks_co)
                    idx_search = close_clicks_co(cl).abs_idx_at_2nd_bat;
                    search_start = (audio_ts(idx_search-search_neighborhood) - co_time_usec)*1e-6;
                    search_end = (audio_ts(idx_search+search_neighborhood) - co_time_usec)*1e-6;
                    search_amplitude = close_clicks_co(cl).snr / mean_snr(i) * mean_snr(3-i);
                    %                 plot([(audio_ts(idx_search)- co_time_usec)*1e-6,(audio_ts(idx_search)- co_time_usec)*1e-6],[-search_amplitude search_amplitude],'k--','Parent',ax{3-i})
                    rectangle('Position',[search_start,-search_amplitude,search_end-search_start,2*search_amplitude],...
                        'FaceColor',[0 0 0 0.2],'EdgeColor','none','Parent',ax{3-i});
                end
                
                co_ind_click = ismember([clicks_struct(i).wrong_other_bat_clicks.peak_ts],audio_ts(co_ind_audio));
                wrong_clicks_co = clicks_struct(i).wrong_other_bat_clicks(co_ind_click);
                plot(([wrong_clicks_co.peak_ts] - co_time_usec)*1e-6,[wrong_clicks_co.snr],'b*','Parent',ax{i});
                
                ylabel('SNR')
                xlabel('time relative to CO (s)')
                legend(p,L{i})
            end
            fig_name = [fig_prefix '_co_dir_' num2str(ii_dir) '_co_' num2str(ii_co) '_self_other_alignment'];
            title(fig_name,'interpreter','none','Parent',ax{1})
            saveas(f,fullfile(fig_folder,fig_name),'png')
        end
    end
    close
    
    
    %% arrange data in matrices
    % bsp
    maxSize = max(cellfun(@numel,bsp_time_to_co));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    %     % a. time to co
    %     rmat = cellfun(fcn,bsp_time_to_co,'UniformOutput',false);
    %     bsp.time_to_co = vertcat(rmat{:});
    % b. relative distance
    rmat = cellfun(fcn,bsp_dis_m_at_co,'UniformOutput',false);
    bsp.dis_m = vertcat(rmat{:});
    
    % clicks
    maxSize = max(cellfun(@numel,clicks_time_to_co));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    %     % a. time to co
    %     % clicks_time_from_co(cellfun(@isempty,clicks_time_from_co)) = {nan};
    %     rmat = cellfun(fcn,clicks_time_to_co,'UniformOutput',false);
    %     clicks.time_to_co = vertcat(rmat{:});
    % b. relative distance
    % clicks_dis_m_at_co(cellfun(@isempty,clicks_dis_m_at_co)) = {nan};
    rmat = cellfun(fcn,clicks_dis_m_at_co,'UniformOutput',false);
    clicks.dis_m = vertcat(rmat{:});
    % c. intensity
    rmat = cellfun(fcn,clicks_intensity_at_co,'UniformOutput',false);
    clicks.intensity = vertcat(rmat{:});
    intensity_vec = clicks.intensity(isfinite(clicks.intensity));
    
    %% calculate calling rate
    % 1. mean:
    %     % a. time to co
    %     bsp_vec_time = bsp.time_to_co(isfinite(bsp.time_to_co));
    %     clicks_vec_time = clicks.time_to_co(isfinite(clicks.time_to_co));
    %
    %     [~, ~, ~, time_to_co_fr, ~,~] ...
    %         = fn_compute_generic_1D_tuning_new_smooth ...
    %         (bsp_vec_time,clicks_vec_time,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
    %
    %     calling_rate.time_to_co = [time_to_co_fr;time_X_bins_vector_of_centers];
    
    % b. relative distance
    bsp_vec_dis = bsp.dis_m(isfinite(bsp.dis_m));
    clicks_vec_dis = clicks.dis_m(isfinite(clicks.dis_m));
    
    [~, ~, ~, dis_m_fr, ~,~] ...
        = fn_compute_generic_1D_tuning_new_smooth ...
        (bsp_vec_dis,clicks_vec_dis,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
    
    calling_rate.dis_m = [dis_m_fr;dis_X_bins_vector_of_centers];
    
    % 2. single CO:
    %     % a. time to co
    %     compute_tuning_curve = @(x,y) fn_compute_generic_1D_tuning_new_smooth ...
    %         (x,y,time_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins/2, frames_per_second, 0,0,0);
    %
    %     [~, ~, ~, single_co_time_to_co_fr, ~,~] = cellfun(compute_tuning_curve,bsp_time_to_co,clicks_time_to_co,'UniformOutput',false);
    %     calling_rate.single_co_time_to_co = cell2mat(single_co_time_to_co_fr);
    
    % b. relative distance
    compute_tuning_curve = @(x,y) fn_compute_generic_1D_tuning_new_smooth ...
        (x,y,dis_X_bins_vector_of_centers, time_spent_minimum_for_1D_bins/2, frames_per_second, 0,0,0);
    
    [~, ~, ~, single_co_dis_m_fr, ~,~] = cellfun(compute_tuning_curve,bsp_dis_m_at_co,clicks_dis_m_at_co,'UniformOutput',false);
    calling_rate.single_co_dis_m = cell2mat(single_co_dis_m_fr);
    
    %% calculate click intensity
    % 1. mean:
    %     % a. time to co
    %     calling_rate.intensity_time_to_co = tuning_curve_by_value_not_rate ...
    %         (clicks_vec_time,bsp_vec_time,time_X_bins_vector,intensity_vec,time_spent_minimum_for_1D_bins,frames_per_second);
    
    % b. relative distance
    calling_rate.intensity_dis_m = tuning_curve_by_value_not_rate ...
        (clicks_vec_dis,bsp_vec_dis,dis_X_bins_vector,intensity_vec,time_spent_minimum_for_1D_bins,frames_per_second);
    
    % 2. single CO:
    %     % a. time to co
    %     fn_tuning_curve_val = @(x,y,z) tuning_curve_by_value_not_rate ...
    %         (x,y,time_X_bins_vector,z,time_spent_minimum_for_1D_bins/2,frames_per_second);
    %
    %     single_co_intensity_time_to_co = cellfun(fn_tuning_curve_val,clicks_time_to_co,bsp_time_to_co,clicks_intensity_at_co,'UniformOutput',false);
    %     calling_rate.single_co_intensity_time_to_co = cell2mat(single_co_intensity_time_to_co);
    
    % b. relative distance
    fn_tuning_curve_val = @(x,y,z) tuning_curve_by_value_not_rate ...
        (x,y,dis_X_bins_vector,z,time_spent_minimum_for_1D_bins/2,frames_per_second);
    
    single_co_intensity_dis_m = cellfun(fn_tuning_curve_val,clicks_dis_m_at_co,bsp_dis_m_at_co,clicks_intensity_at_co,'UniformOutput',false);
    calling_rate.single_co_intensity_dis_m = cell2mat(single_co_intensity_dis_m);
    
    
    %% save struct
    co(ii_dir).clicks = clicks;
    co(ii_dir).bsp = bsp;
    co(ii_dir).calling_rate = calling_rate;
    
end
end