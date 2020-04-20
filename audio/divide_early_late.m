function calling_rate = divide_early_late(clicks_struct_name,click_rate_folder,name_prefix)

load(clicks_struct_name,'clicks_struct')

diff_th=1;
smooth_wind = 5;
smooth_wind_diff = 9;
cr_upsample = 100;
solo_prcnt = [80 20];

for ii_dir = 1:2
    %             mean_cr = clicks_struct(1).co(ii_dir).calling_rate.dis_m(1,:);
    %             mean_cr_diff = diff(mean_cr);
    x_bins = clicks_struct(1).co(ii_dir).calling_rate.dis_m(2,:);
    single_co_cr = clicks_struct(1).co(ii_dir).calling_rate.single_co_dis_m;
    single_co_smooth = smoothdata(single_co_cr,2,'movmean',smooth_wind);
    %             single_co_smooth_for_diff = smoothdata(single_co_cr,2,'movmean',smooth_wind_diff);
    
    %upsampling the calling rate
    x_bins_up = linspace(x_bins(1),x_bins(end),length(x_bins)*cr_upsample);
    single_co_smooth_up = (interp1(x_bins',single_co_smooth',x_bins_up'))';
    
    %             all_cr = [mean_cr;single_co_cr];
    %             all_cr_diff = [mean_cr_diff;single_co_cr_diff];
    [max_rate,max_rate_ind] = max(single_co_smooth_up,[],2);
    cr_bl = single_co_smooth_up(:,1);
    %             peak_height = max_rate - cr_bl;
    high_bl = cr_bl >= 0.5*max_rate;
    
    nbins = length(x_bins_up);
    before_peak = cellfun(@(x) [true(1,x),false(1,nbins-x)],num2cell(max_rate_ind),'UniformOutput',false);
    %             single_co_smooth(high_bl,:) = [];
    %             single_co_smooth_for_diff(high_bl,:) = [];
    %             single_co_cr_diff = diff(single_co_smooth_for_diff,1,2);
    %             rate_inc_deriv_cell = cellfun(@(x) find(x>=diff_th,1),num2cell(single_co_cr_diff,2),'UniformOutput',false);
    %             rate_inc_half_height_cell = cellfun(@(x,y) find(x>=0.5*y,1),num2cell(single_co_smooth,2),num2cell(max_rate(~high_bl)),'UniformOutput',false);
    
    rate_inc_half_height_cell = cellfun(@(x,y) x<0.5*y,num2cell(single_co_smooth_up,2),num2cell(max_rate),'UniformOutput',false);
    rate_inc_half_height_cell = cellfun(@(x,y) find(x & y,1,'last'),rate_inc_half_height_cell,before_peak,'UniformOutput',false);
    %             problem_with_deriv = cellfun(@isempty,rate_inc_deriv_cell);
    problem_with_half_height = cellfun(@isempty,rate_inc_half_height_cell);
    %             rate_inc_deriv_cell(problem_with_deriv) = num2cell(NaN);
    rate_inc_half_height_cell(problem_with_half_height) = num2cell(NaN);
    %             rate_inc_deriv = cell2mat(rate_inc_deriv_cell);
    rate_inc_half_height = cell2mat(rate_inc_half_height_cell);
    %             thresh_deriv = floor(nanmedian(rate_inc_deriv));
    thresh_half_height = floor(nanmedian(rate_inc_half_height));
    
    solo_cr = clicks_struct(1).solo(ii_dir).calling_rate;
    mean_solo_cr = mean(solo_cr);
    ci_solo_cr = [prctile(solo_cr,solo_prcnt(1)),prctile(solo_cr,solo_prcnt(2))];
    err = [ci_solo_cr(1)-mean_solo_cr;mean_solo_cr-ci_solo_cr(2)];
    
    rate_inc_solo_cell = cellfun(@(x) x<ci_solo_cr(1),num2cell(single_co_smooth_up,2),'UniformOutput',false);
    rate_inc_solo_cell = cellfun(@(x,y) find(x & y,1,'last'),rate_inc_solo_cell,before_peak,'UniformOutput',false);
    problem_with_solo = cellfun(@isempty,rate_inc_solo_cell);
    rate_inc_solo_cell(problem_with_solo) = num2cell(NaN);
    rate_inc_solo = cell2mat(rate_inc_solo_cell);
    thresh_solo = floor(nanmedian(rate_inc_solo));
    
    %% single co figure
    click_dis = clicks_struct(1).co(ii_dir).clicks.dis_m;
    %             click_dis(high_bl,:) = [];
    figure
    for ico=1:length(max_rate)
        if high_bl(ico)
            continue
        end
        clf
        axes('Position',[0.1 0.2 0.8 0.7])
        hold on
        xlabel('inter-bat distance (m)')
        %                 yyaxis left
        shadedErrorBar([x_bins(1) x_bins(end)],repmat(mean_solo_cr,1,2),repmat(err,1,2));
        plot(x_bins,single_co_smooth(ico,:),'LineWidth',2,'color',[0 0.4470 0.7410])
        xlim([x_bins(1) x_bins(end)])
        ylim([0 35])
        ylabel('Calling Rate (Hz)')
        %                 yyaxis right
        %                 plot(x_bins(2:end),single_co_cr_diff(ico,:),'LineWidth',2,'color',[0.8500 0.3250 0.0980])
        %                 ylim([-7 7])
        %                 if ~problem_with_deriv(ico)
        %                     xline(x_bins(rate_inc_deriv(ico)),'--','color',[0.8500 0.3250 0.0980]);
        %                 end
        if ~problem_with_half_height(ico)
            xline(x_bins_up(rate_inc_half_height(ico)),'--','color',[0 0.4470 0.7410]);
        end
        if ~problem_with_solo(ico)
            xline(x_bins_up(rate_inc_solo(ico)),'--','color',[0.5 0.5 0.5]);
        end
        if problem_with_solo(ico) && problem_with_half_height(ico)
            text(2,5,'problem with both threshold methods')
        end
        %                 ylabel('Calling Rate diff (Hz)')
        title({name_prefix,['dir_' num2str(ii_dir) '_co_' num2str(ico)]},'Interpreter','none')
        legend(['Solo calling rate (Hz) - mean+' num2str(solo_prcnt(2)) '-' num2str(solo_prcnt(1)) '%'])
        
        axes('Position',[0.1 0.05 0.8 0.05])
        plot([click_dis(ico,:);click_dis(ico,:)],repmat([1;-1],1,length(click_dis(ico,:))),'b');
        axis off
        xlim([x_bins(1) x_bins(end)])
        ylim([-3 3])
        
        saveas(gcf,fullfile(click_rate_folder,'single_co',[name_prefix '_dir_' num2str(ii_dir) '_co_' num2str(ico)]),'png')
    end
    close
    
    %% all COs figure
    figure('WindowState','maximized')
    for i=0:1
        clf
        if i
            thresh_used = thresh_half_height;
            co_early = rate_inc_half_height < thresh_used;
            rate_inc = rate_inc_half_height;
            prob = problem_with_half_height;
            criterion_used = 'half-height criterion';
            
            %                     thresh_used = thresh_deriv;
            %                     co_early = rate_inc_deriv < thresh_used;
            %                     rate_inc = rate_inc_deriv;
            %                     prob = problem_with_deriv;
            %                     criterion_used = ['diff criterion; diff thresh = ' num2str(diff_th)];
        else
            thresh_used = thresh_solo;
            co_early = rate_inc_solo < thresh_used;
            rate_inc = rate_inc_solo;
            prob = problem_with_solo;
            criterion_used = 'solo criterion';
        end
        hold on
        plot(x_bins,single_co_smooth(co_early & ~prob,:),'color',[1 0 0 .5],'LineWidth',0.7)
        plot(x_bins,single_co_smooth(~co_early & ~prob,:),'color',[0 0 1 .5],'LineWidth',0.7)
        plot(x_bins,nanmean(single_co_smooth(co_early & ~prob,:)),'color',[1 0 0 1],'LineWidth',2.5)
        plot(x_bins,nanmean(single_co_smooth(~co_early & ~prob,:)),'color',[0 0 1 1],'LineWidth',2.5)
        plot(x_bins(rate_inc(co_early & ~prob)),cellfun(@(x,y) x(y),num2cell(single_co_smooth(co_early & ~prob,:),2),num2cell(rate_inc(co_early & ~prob))),'ro')
        plot(x_bins(rate_inc(~co_early & ~prob)),cellfun(@(x,y) x(y),num2cell(single_co_smooth(~co_early & ~prob,:),2),num2cell(rate_inc(~co_early & ~prob))),'bo')
        xline(x_bins(thresh_used),'--');
        xlabel('inter-bat distance (m)')
        ylabel('Calling Rate (Hz)')
        title({name_prefix,'CO events divided by calling-rate increase latency',['dir_' num2str(ii_dir) '; ' criterion_used]},'Interpreter','none')
        axes('Position',[.7 .7 .2 .2])
        hold on
        cla
        [~,edges] = histcounts(x_bins_up(rate_inc));
        histogram(x_bins_up(rate_inc(co_early & ~prob)),edges,'FaceColor',[1 0 0])
        histogram(x_bins_up(rate_inc(~co_early & ~prob)),edges,'FaceColor',[0 0 1])
        xlabel('inter-bat distance (m)')
        saveas(gcf,fullfile(click_rate_folder,['cr_hh_' num2str(i) '_' name_prefix,'_dir_' num2str(ii_dir) '.png']))
    end
    close
    
    %% click amplitude figure
    solo_amp = clicks_struct(1).solo(ii_dir).click_intensity;
    mean_solo_amp = nanmean(solo_amp);
    err = [prctile(solo_amp,90)-mean_solo_amp;mean_solo_amp-prctile(solo_amp,10)];
    
    co_amp = clicks_struct(1).co(ii_dir).calling_rate.single_co_intensity_dis_m;
    co_amp_mean = clicks_struct(1).co(ii_dir).calling_rate.intensity_dis_m;
    co_amp_sm=smoothdata(co_amp,2,'movmean',smooth_wind);
    co_amp_sm(isnan(co_amp)) = NaN;
    co_amp_mean=smoothdata(co_amp_mean,'movmean',smooth_wind);
    figure('WindowState','maximized')
    hold on
    p_solo=shadedErrorBar([x_bins(1) x_bins(end)],repmat(mean_solo_amp,1,2),repmat(err,1,2));
    set(get(get(p_solo.mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    p_co=plot(x_bins,co_amp,'color',[0.3 0 0.3 .5],'LineWidth',0.7,'HandleVisibility','off');
    plot(x_bins,co_amp_mean,'color',[0.3 0 0.3 1],'LineWidth',2.5)
    xlabel('inter-bat distance (m)')
    ylabel('SNR')
    title({name_prefix,'Click Amplitude',['dir_' num2str(ii_dir)]},'Interpreter','none')
    legend({'Solo','CO'})
    saveas(gcf,fullfile(click_rate_folder,['amp_' name_prefix,'_dir_' num2str(ii_dir) '.png']))
    close
    
    %% save co_early and other params in the struct!!! (for the 2 directions)
    %calling_rate(ii_dir) = ;
end