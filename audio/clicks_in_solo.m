function solo = clicks_in_solo(bsp_proc_data,behavioral_modes,tag_i,solo_param_file_name,clicks_struct,plot_flight_flag,audio_filt_struct,fig_folder,fig_prefix,th)
load(solo_param_file_name)
us_factor=1e6;
%% create vectors of variables

% a. time stamps
clicks_ts = [clicks_struct(strcmp({clicks_struct.bat},'self')).clusters.peak_ts];
if plot_flight_flag
audio_ts = audio_filt_struct(1).ts;
end
% intensity
click_intensity = [clicks_struct(strcmp({clicks_struct.bat},'self')).clusters.snr];

solo_ind = behavioral_modes.solo_ind;

%% loop for the 2 directions
for ii_dir = 1:2
    
    dir_ind = behavioral_modes.directional_ind{ii_dir};    
    %find bsp parameters
    bsp_during_solo_ts_usec = bsp_proc_data(tag_i).ts(intersect(solo_ind,dir_ind));
    bsp_during_solo_x_pos = bsp_proc_data(tag_i).pos(intersect(solo_ind,dir_ind),1);
    
    % create a matrix from the data, such as each row is one flight
    new_flight_dis_criteria = dis_criteria(ii_dir);
    new_flight_dis_ind = abs(diff(bsp_during_solo_x_pos)) > new_flight_dis_criteria;
    new_flight_time_ind = diff(bsp_during_solo_ts_usec) > new_flight_time_criteria;
    new_flight_ind = find(new_flight_dis_ind | new_flight_time_ind);
    new_flight_start_ind = [1;new_flight_ind+1];
    new_flight_end_ind = [new_flight_ind;length(bsp_during_solo_ts_usec)];
    n_flights = length(new_flight_start_ind);
    
    % prepare empty variables
    bsp_ts_usec = cell(n_flights,1);
    bsp_x_pos = cell(n_flights,1);
    clicks_ts_usec = cell(n_flights,1);
    clicks_x_pos = cell(n_flights,1);
    calling_rate = zeros(n_flights,1);
    click_intensity_flight = zeros(n_flights,1);
    
    calculate_rate = @(x,y) length(x) ./ (y(end) - y(1)) * us_factor;
    
    % loop for every flight
    f=figure('WindowState','maximized');
    for ii_flight = 1:n_flights
        flight_ind = new_flight_start_ind(ii_flight):new_flight_end_ind(ii_flight);
        bsp_ts_usec{ii_flight} = bsp_during_solo_ts_usec(flight_ind)';
        bsp_x_pos{ii_flight} = bsp_during_solo_x_pos(flight_ind)';
        clicks_ind = clicks_ts >= bsp_ts_usec{ii_flight}(1) & clicks_ts <= bsp_ts_usec{ii_flight}(end);
        clicks_ts_usec{ii_flight} = clicks_ts(clicks_ind);
        clicks_x_pos{ii_flight} = interp1(bsp_ts_usec{ii_flight},bsp_x_pos{ii_flight},clicks_ts_usec{ii_flight});
        
        calling_rate(ii_flight) = calculate_rate(clicks_ts_usec{ii_flight},bsp_ts_usec{ii_flight});
        click_intensity_flight(ii_flight) = mean(click_intensity(clicks_ind));
        
        if plot_flight_flag
        % plot audio signal for this flight
        clf
        flight_ind_audio = audio_ts >= bsp_ts_usec{ii_flight}(1) & audio_ts <= bsp_ts_usec{ii_flight}(end);
        audio_ts_flight = audio_ts(flight_ind_audio);
        ts_0 = audio_ts_flight(1);
        C = {'g','r'};
        L = {'self bat logger (implanted)','other bat logger (not implanted)'};
        i=1;
        hold on
        flight_ind_click = ismember([clicks_struct(i).clusters.peak_ts],audio_ts_flight);
        plot((audio_ts_flight - ts_0)*1e-6,audio_filt_struct(i).filt(flight_ind_audio)./clicks_struct(i).aud_BL,C{i})
        p=plot(([clicks_struct(i).clusters(flight_ind_click).peak_ts] - ts_0)*1e-6,[clicks_struct(i).clusters(flight_ind_click).snr],'k*');
        yline(th,'--');
        yline(-th,'--');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        ylabel('SNR')
        xlabel('time in flight (s)')
        fig_name = [fig_prefix '_solo_dir_' num2str(ii_dir) '_flight_' num2str(ii_flight)];
        title({fig_name,L{i},['solo calling rate: ' num2str(calling_rate(ii_flight)) ' Hz']},'interpreter','none')
        saveas(f,fullfile(fig_folder,fig_name),'png')
        end
    end
    close
    
    %% calculate calling rate in each flight
%     calling_rate = cellfun(@(x,y) length(x) ./ (y(end) - y(1)) * us_factor,clicks_ts_usec,bsp_ts_usec);
     
%     solo(ii_dir).clicks = clicks;
%     solo(ii_dir).bsp = bsp;
    solo(ii_dir).calling_rate = calling_rate;
    solo(ii_dir).click_intensity = click_intensity_flight;
    
end
