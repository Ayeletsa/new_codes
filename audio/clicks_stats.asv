% plot click statistics

if ~exist('fig_prefix','var')
    fig_prefix = clicks_struct(1).name;
end

%% plot click statistics
figure('WindowState','maximized')
C = {'g','r'};
sgtitle(fig_prefix,'interpreter','none')
for ibat=1:2
    ibat_idx = 3*(ibat-1)+1;
    
    subplot(4,6,ibat_idx+0)
    clicks_snr = [clicks_struct(ibat).clusters.snr];
    histogram(clicks_snr,'FaceColor',C{ibat});
    title([clicks_struct(ibat).bat ' - snr'])
    
    subplot(4,6,ibat_idx+1)
    clicks_length = [clicks_struct(ibat).clusters.Area]*1e3/fs_aud;
    histogram(clicks_length,'FaceColor',C{ibat});
    xlabel('ms')
    title([clicks_struct(ibat).bat ' - length'])
    
    subplot(4,6,ibat_idx+2)
    idx = [clicks_struct(ibat).clusters.peak_abs_idx] - [clicks_struct(ibat).clusters.start] + 1;
    histogram(idx*1e3/fs_aud,'FaceColor',C{ibat});
    xlabel('ms')
    title([clicks_struct(ibat).bat ' - rise time'])
    
    subplot(4,6,ibat_idx+6)
    scatter(clicks_length,clicks_snr,5,C{ibat},'marker','.')
    xlabel('length(ms)')
    ylabel('snr')
    title([clicks_struct(ibat).bat ' - snr vs. length'])
    
    subplot(4,6,ibat_idx+12)
    mid_low_spect_ratio = [clicks_struct(ibat).clusters.spectrum_ratio_mid_low];
    histogram(mid_low_spect_ratio,'binWidth',1,'FaceColor',C{ibat});
    xlabel('ratio (db)')
    title([clicks_struct(ibat).bat ' - mid/low spectral ratio'])
    
    subplot(4,6,ibat_idx+13)
    scatter(mid_low_spect_ratio,clicks_snr,5,C{ibat},'marker','.')
    xlabel('ratio (db)')
    ylabel('snr')
    title([clicks_struct(ibat).bat ' - snr vs. mid/low spectral ratio'])
    
    subplot(4,6,ibat_idx+14)
    scatter(mid_low_spect_ratio,clicks_length,5,C{ibat},'marker','.')
    xlabel('ratio (db)')
    ylabel('length (ms)')
    title([clicks_struct(ibat).bat ' - length vs. mid/low spectral ratio'])
    
    subplot(4,6,ibat_idx+18)
%     hold on
%     high_mid_spect_ratio = [clicks_struct(ibat).clusters.spectrum_ratio_high_mid];
%     high_mid_spect_ratio_noise = [clicks_struct(ibat).wrong_mid_low_spectrum.spectrum_ratio_high_mid];
%     histogram(high_mid_spect_ratio,'binWidth',1,'FaceColor',C{ibat},'faceAlpha',0.6)
%     histogram(high_mid_spect_ratio_noise,'binWidth',1,'FaceColor',[0.5 0.5 0.5])
%     title([clicks_struct(ibat).bat ' - high/mid spectral ratio'])
    
    ici = diff([clicks_struct(ibat).clusters.peak_abs_idx])/us2fs/1e3;
    histogram(ici(ici<50),'binWidth',1,'FaceColor',C{ibat})
    xlabel('ms')
    title([clicks_struct(ibat).bat ' - ICI'])
    
end
saveas(gcf,fullfile(fig_folder,[fig_prefix '_click_stats.png']))
close

