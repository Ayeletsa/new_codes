function calling_rate = divide_early_late(clicks_struct_name,click_rate_folder,name_prefix,co_param_file_name,audio_param_file_name,early_late_file_name)

load(clicks_struct_name,'clicks_struct')
load(co_param_file_name)
load(audio_param_file_name)
load(early_late_file_name)

name_prefix = strrep(name_prefix,'.mat','');

% smooth_wind = 5;
% cr_upsample = 100;
% solo_prcnt = [75 25];
% count_dc_once = false;

x_bins = dis_X_bins_vector_of_centers;
x_bins_up = linspace(x_bins(1),x_bins(end),length(x_bins)*cr_upsample);
nbins = length(x_bins_up);

x_amp_bins = amp_dis_X_bins_vector_of_centers;

for ii_dir = 1:2
    if count_dc_once
        co_cr = clicks_struct(1).co(ii_dir).calling_rate.single_co_dis_m_once;
        solo_cr = clicks_struct(1).solo(ii_dir).calling_rate_once;
        
    else    
        co_cr = clicks_struct(1).co(ii_dir).calling_rate.single_co_dis_m;
        solo_cr = clicks_struct(1).solo(ii_dir).calling_rate;
    end
    
    %smoothing and upsampling the calling rate
    co_smooth = smoothdata(co_cr,2,'movmean',smooth_wind);
    co_smooth_up = (interp1(x_bins',co_smooth',x_bins_up'))';
    
    single_co_up = num2cell(co_smooth_up,2);
    
    %calculate when the rate increase
    [max_rate,max_rate_ind] = max(co_smooth_up,[],2);
    bl_ind = cellfun(@(x) find(~isnan(x),1),single_co_up);
    cr_bl = cellfun(@(x,y) x(y),single_co_up,num2cell(bl_ind));
    high_bl = cr_bl >= 0.5*max_rate;
    
    before_peak = cellfun(@(x) [true(1,x),false(1,nbins-x)],num2cell(max_rate_ind),'UniformOutput',false);
    
    % half-height criterion
    rate_inc_half_height_cell = cellfun(@(x,y) x<0.5*y,single_co_up,num2cell(max_rate),'UniformOutput',false);
    rate_inc_half_height_cell = cellfun(@(x,y) find(x & y,1,'last'),rate_inc_half_height_cell,before_peak,'UniformOutput',false);
    rate_inc_half_height_cell(high_bl) = num2cell(NaN);
    rate_inc_half_height = cell2mat(rate_inc_half_height_cell);
    thresh_half_height = floor(nanmedian(rate_inc_half_height));
    
    % solo criterion
    mean_solo_cr = nanmean(solo_cr);
    ci_solo_cr = [prctile(solo_cr,solo_prcnt(1)),prctile(solo_cr,solo_prcnt(2))];
    err_solo = [ci_solo_cr(1)-mean_solo_cr;mean_solo_cr-ci_solo_cr(2)];
    
    rate_inc_solo_cell = cellfun(@(x) x<ci_solo_cr(1),single_co_up,'UniformOutput',false);
    rate_inc_solo_cell = cellfun(@(x,y) find(x & y,1,'last'),rate_inc_solo_cell,before_peak,'UniformOutput',false);
    high_attention_by_solo = cellfun(@isempty,rate_inc_solo_cell);
    low_attention_by_solo = max_rate < ci_solo_cr(1);
    rate_inc_solo_cell(high_attention_by_solo | low_attention_by_solo | high_bl) = num2cell(NaN);
    rate_inc_solo = cell2mat(rate_inc_solo_cell);
    thresh_solo = floor(nanmedian(rate_inc_solo));
    
    calling_rate(ii_dir).th_half_height = thresh_half_height;
    calling_rate(ii_dir).half_height_early = rate_inc_half_height < thresh_half_height & ~high_bl;    
    calling_rate(ii_dir).half_height_late = rate_inc_half_height >= thresh_half_height & ~high_bl;
    calling_rate(ii_dir).th_solo = thresh_solo;
    calling_rate(ii_dir).solo_early = rate_inc_solo < thresh_solo & ~high_attention_by_solo & ~low_attention_by_solo & ~high_bl;    
    calling_rate(ii_dir).solo_late = rate_inc_solo >= thresh_solo & ~high_attention_by_solo & ~low_attention_by_solo & ~high_bl;
    calling_rate(ii_dir).high_bl = high_bl;
    calling_rate(ii_dir).high_attention_by_solo = high_attention_by_solo;
    calling_rate(ii_dir).low_attention_by_solo = low_attention_by_solo;
    
    
%     %% single co figure
%     click_dis = clicks_struct(1).co(ii_dir).clicks.dis_m;
%     single_clicks = clicks_struct(1).co(ii_dir).clicks.not_DC;
%     click_dis_single = click_dis .* single_clicks;
%     
%     figure
%     C = [0 0.4470 0.7410];%,[0.8500 0.3250 0.0980]};
%     for ico=1:length(max_rate)
%         if high_bl(ico)
%             continue
%         end
%         clf
%         
%         axes('Position',[0.1 0.2 0.8 0.7])
%         hold on
%         xlabel('inter-bat distance (m)')
%             p=shadedErrorBar([x_bins(1) x_bins(end)],repmat(mean_solo_cr,1,2),repmat(err_solo,1,2));
%             set(p.edge,'visible','off')
%             p.patch.FaceColor = C;
%             plot(x_bins,co_smooth(ico,:),'LineWidth',2,'color',C)
% 
%             if ~isnan(rate_inc_solo(ico))
%                 xline(x_bins_up(rate_inc_solo(ico)),'--','color',[0.5 0.5 0.5]);
%             elseif high_attention_by_solo(ico)
%                 text(2,5,'high attention by solo')
%             elseif low_attention_by_solo(ico)
%                 text(2,5,'low attention by solo') 
%             end
% 
%         xlim([x_bins(1) x_bins(end)])
%         ylim([0 35])
%         ylabel('Calling Rate (Hz)')
%         title({name_prefix,['dir_' num2str(ii_dir) '_co_' num2str(ico)]},'Interpreter','none')
%         legend(['Solo calling rate (Hz) - mean+' num2str(solo_prcnt(2)) '-' num2str(solo_prcnt(1)) '%'])
%         
%         axes('Position',[0.1 0.05 0.8 0.05])
%         hold on
%         plot([click_dis(ico,:);click_dis(ico,:)],repmat([1;-1],1,length(click_dis(ico,:))),'b');
%         plot([click_dis_single(ico,:);click_dis_single(ico,:)],repmat([1;-1],1,length(click_dis_single(ico,:))),'r');
%         axis off
%         xlim([x_bins(1) x_bins(end)])
%         ylim([-3 3])
%         
%         saveas(gcf,fullfile(click_rate_folder,'single_co',[name_prefix '_dir_' num2str(ii_dir) '_co_' num2str(ico)]),'png')
%     end
%     close
    
    %% all COs figure
    figure('WindowState','maximized')
    for i=0:1
        clf
        if i
            thresh_used = thresh_half_height;
            co_early = calling_rate(ii_dir).half_height_early;
            co_late = calling_rate(ii_dir).half_height_late;
            rate_inc = rate_inc_half_height;
            
            criterion_used = 'half-height criterion';
        else
            thresh_used = thresh_solo;
            co_early = calling_rate(ii_dir).solo_early;
            co_late = calling_rate(ii_dir).solo_late;
            rate_inc = rate_inc_solo;
            criterion_used = 'solo criterion';
        end
        
        hold on
        plot(x_bins,co_smooth(co_early,:),'color',[1 0 0 .5],'LineWidth',0.7)
        plot(x_bins,co_smooth(co_late,:),'color',[0 0 1 .5],'LineWidth',0.7)
        plot(x_bins,nanmean(co_smooth(co_early,:)),'color',[1 0 0 1],'LineWidth',2.5)
        plot(x_bins,nanmean(co_smooth(co_late,:)),'color',[0 0 1 1],'LineWidth',2.5)
        
        rate_at_inc = nan(length(rate_inc),1);
        for ico=1:length(rate_inc)
            if ~isnan(rate_inc(ico))
                rate_at_inc(ico) = co_smooth_up(ico,rate_inc(ico));
            end
        end
        plot(x_bins_up(rate_inc(co_early)),rate_at_inc(co_early),'ro')
        plot(x_bins_up(rate_inc(co_late)),rate_at_inc(co_late),'bo')
        
        xline(x_bins_up(thresh_used),'--');
        xlabel('inter-bat distance (m)')
        ylabel('Calling Rate (Hz)')
        title({name_prefix,'CO events divided by calling-rate increase latency',['dir_' num2str(ii_dir) '; ' criterion_used]},'Interpreter','none')
        
        axes('Position',[.7 .7 .2 .2])
        hold on
        cla
        histogram(x_bins_up(rate_inc(co_early)),'BinWidth',2,'FaceColor',[1 0 0])
        histogram(x_bins_up(rate_inc(co_late)),'BinWidth',2,'FaceColor',[0 0 1])
        xlabel('inter-bat distance (m)')
        saveas(gcf,fullfile(click_rate_folder,['cr_hh_' num2str(i) '_' name_prefix,'_dir_' num2str(ii_dir) '.png']))
    end
    close
    
%     %% click amplitude figure
%     solo_amp = clicks_struct(1).solo(ii_dir).click_intensity;
%     median_solo_amp = nanmedian(solo_amp);
%     err_solo_amp = [prctile(solo_amp,solo_prcnt(1))-median_solo_amp;median_solo_amp-prctile(solo_amp,solo_prcnt(2))];
%     
%     co_amp = clicks_struct(1).co(ii_dir).calling_rate.single_co_intensity_dis_m;
%     co_amp_median = nanmedian(co_amp,1);
%     err_co = [prctile(co_amp,solo_prcnt(1),1)-co_amp_median;co_amp_median-prctile(co_amp,solo_prcnt(2),1)];
% %     co_amp_mean = clicks_struct(1).co(ii_dir).calling_rate.intensity_dis_m;
% %     co_amp_sm=smoothdata(co_amp,2,'movmean',smooth_wind);
% %     co_amp_sm(isnan(co_amp)) = NaN;
%     co_amp_median=smoothdata(co_amp_median,'movmean',smooth_wind);
%     figure('WindowState','maximized')
%     hold on
%     p_solo=shadedErrorBar([x_bins(1) x_bins(end)],repmat(median_solo_amp,1,2),repmat(err_solo_amp,1,2));
%     set(get(get(p_solo.mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% %     p_co=plot(x_bins,co_amp_sm,'color',[0.3 0 0.3 .5],'LineWidth',0.7,'HandleVisibility','off');
% %     plot(x_bins,co_amp_mean,'color',[0.3 0 0.3 1],'LineWidth',2.5)
%     p_solo.mainLine.LineWidth = 2.5;
%     set(p_solo.edge,'visible','off')
%     p_co=shadedErrorBar(x_amp_bins,co_amp_median,err_co);
%     p_co.patch.FaceColor = [0.3 0 0.3];
%     set(p_co.mainLine,'LineWidth',2.5,'Color',[0.3 0 0.3])
%     set(p_co.edge,'visible','off')
%     xlabel('inter-bat distance (m)')
%     ylabel('SNR')
%     title({name_prefix,'Click Amplitude',['dir_' num2str(ii_dir)]},'Interpreter','none')
%     L=legend([p_solo.mainLine,p_co.mainLine],{'Solo','CO'});
%     title(L,['median + ' num2str(solo_prcnt(2)) '-' num2str(solo_prcnt(1)) '%'])
%     saveas(gcf,fullfile(click_rate_folder,['amp_' name_prefix,'_dir_' num2str(ii_dir) '.png']))
%     close
    
end