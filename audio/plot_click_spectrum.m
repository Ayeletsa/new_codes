% plot weird clicks with their spectral properties
figure('WindowState','maximized')
max_per_fig = 5;
margin_close = 3;%ms
margin_far = 300;%ms

min_snr = 0;
max_ratio = -17;

for ibat = 1%:2
clicks_to_use = clicks_struct(ibat).clusters([clicks_struct(ibat).clusters.snr] > min_snr & [clicks_struct(ibat).clusters.spectrum_ratio_mid_low] < max_ratio);
ncl = length(clicks_to_use);

q = floor(ncl/max_per_fig);
r = rem(ncl,max_per_fig);
v = [repmat(max_per_fig,1,q),r];
ve = cumsum(v);
vs = ve - v + 1;

fig_name = [fig_prefix '_' clicks(ibat).bat '_largeSNR_low_spect'];

for ifig=1:length(v)
    if v(ifig)==0
        break
    end
    clf
    cl_vec = vs(ifig):ve(ifig);
    for icl=1:length(cl_vec)
        cl = cl_vec(icl);
        
        % zoom out on the filtered signal
        subplot(4,max_per_fig,icl)
        plot_wind_far = -1e3*margin_far*us2fs:1e3*margin_far*us2fs;
        ind_to_plot_far = clicks_to_use(cl).peak_abs_idx + plot_wind_far;
        plot(plot_wind_far*1e3/fs_aud,clicks(ibat).filt(ind_to_plot_far)/clicks(ibat).aud_BL)
        xlabel('ms')
        ylabel('snr')
        ylim([-2000 2000])
        xlim([-margin_far margin_far])
        title(num2str(cl))
        
        % close look on the unfiltered signal
        subplot(4,max_per_fig,icl+max_per_fig)
        ind_fft = clicks_to_use(cl).peak_abs_idx + fft_window;
        plot_wind_close = -1e3*margin_close*us2fs:1e3*margin_close*us2fs;
        ind_to_plot = clicks_to_use(cl).peak_abs_idx + plot_wind_close;
        plot(plot_wind_close*1e3/fs_aud,clicks(ibat).unfilt(ind_to_plot)/clicks(ibat).aud_BL)
        xline(fft_window(1)*1e3/fs_aud,'r--');
        xline(fft_window(end)*1e3/fs_aud,'r--');
        ylim([-2000 2000])
        xlabel('ms')
        ylabel('snr')
        xlim([-margin_close margin_close])
        title(['snr=' num2str(clicks_to_use(cl).snr)])
        
        
        % fft spectrum
        subplot(4,max_per_fig,icl+2*max_per_fig)
        fft_click = fft(clicks(ibat).unfilt(ind_fft));
        fft_click = abs(fft_click/L);
        fft_click_single_side = fft_click(1:L/2+1);
        fft_click_single_side(2:end-1) = 2*fft_click_single_side(2:end-1);
        f_fft = (fs_aud*(0:(L/2))/L)*10^-3; %KHz;
        plot(f_fft,fft_click_single_side)
        xline(low_freq_band(1),'g--');
        xline(low_freq_band(2),'g--');
        xline(high_freq_band(1),'m--');
        xline(high_freq_band(2),'m--');
        xlim([0 50])
        xlabel('kHz')
        ylabel('magnitude')
        title(['ratio=' num2str(clicks_to_use(cl).spectrum_ratio_mid_low) 'dB'])
        
        % spectrogram
        subplot(4,max_per_fig,icl+3*max_per_fig)
        [s_stft,f_stft,t_stft] = stft(clicks(ibat).unfilt(ind_to_plot),fs_aud);
        s_stft = mag2db(abs(s_stft));
        c_limits = [max(s_stft(:))-60,max(s_stft(:))];
        imagesc((t_stft*1e3)+plot_wind_close(1)*1e3/fs_aud,f_stft*1e-3,s_stft,c_limits)
        ylim([0 50])
        axis xy
        xline(fft_window(1)*1e3/fs_aud,'r--');
        xline(fft_window(end)*1e3/fs_aud,'r--');
        yline(low_freq_band(1),'g--');
        yline(low_freq_band(2),'g--');
        yline(mid_freq_band(1),'m--');
        yline(mid_freq_band(2),'m--');
        c = colorbar;
        c.Label.String = 'db';
        xlabel('ms')
        ylabel('kHz')
    end
    sgtitle(fig_name,'interpreter','none')
    saveas(gcf,fullfile(fig_folder,[fig_name '_' num2str(ifig) '.png']))
end
end
close
 