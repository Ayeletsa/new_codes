function calling_rate = divide_early_late(clicks_struct_name,click_rate_folder,name_prefix,co_param_file_name)

load(clicks_struct_name,'clicks_struct')
load(co_param_file_name)

name_prefix = strrep(name_prefix,'.mat','');

diff_th=1;
smooth_wind = 5;
smooth_wind_diff = 9;
cr_upsample = 100;
solo_prcnt = [75 25];

for ii_dir = 1:2
    for jj=1:2
        
    %             mean_cr = clicks_struct(1).co(ii_dir).calling_rate.dis_m(1,:);
    %             mean_cr_diff = diff(mean_cr);
    x_bins = dis_X_bins_vector_of_centers;
    
    if jj==1
        single_co_cr = clicks_struct(1).co(ii_dir).calling_rate.single_co_dis_m;
        solo_cr = clicks_struct(1).solo(ii_dir).calling_rate;
    else
        single_co_cr = clicks_struct(1).co(ii_dir).calling_rate.single_co_dis_m_once;
        solo_cr = clicks_struct(1).solo(ii_dir).calling_rate_once;
    end
    single_co_smooth{jj} = smoothdata(single_co_cr,2,'movmean',smooth_wind);
    %             single_co_smooth_for_diff = smoothdata(single_co_cr,2,'movmean',smooth_wind_diff);
    
    %upsampling the calling rate
    x_bins_up = linspace(x_bins(1),x_bins(end),length(x_bins)*cr_upsample);
    single_co_smooth_up = (interp1(x_bins',single_co_smooth{jj}',x_bins_up'))';
    
    %             all_cr = [mean_cr;single_co_cr];
    %             all_cr_diff = [mean_cr_diff;single_co_cr_diff];
    [max_rate,max_rate_ind] = max(single_co_smooth_up,[],2);
    cr_bl = cellfun(@(x) x(find(~isnan(x),1)),num2cell(single_co_smooth_up,2));
%     cr_bl = single_co_smooth_up(:,1);
    %             peak_height = max_rate - cr_bl;
    high_bl{jj} = cr_bl >= 0.5*max_rate;
    
    nbins = length(x_bins_up);
    before_peak = cellfun(@(x) [true(1,x),false(1,nbins-x)],num2cell(max_rate_ind),'UniformOutput',false);
    %             single_co_smooth(high_bl,:) = [];
    %             single_co_smooth_for_diff(high_bl,:) = [];
    %             single_co_cr_diff = diff(single_co_smooth_for_diff,1,2);
    %             rate_inc_deriv_cell = cellfun(@(x) find(x>=diff_th,1),num2cell(single_co_cr_diff,2),'UniformOutput',false);
    %             rate_inc_half_height_cell = cellfun(@(x,y) find(x>=0.5*y,1),num2cell(single_co_smooth,2),num2cell(max_rate(~high_bl)),'UniformOutput',false);
    
%     problem_with_deriv = cellfun(@isempty,rate_inc_deriv_cell);
%     rate_inc_deriv_cell(problem_with_deriv) = num2cell(NaN);
%     rate_inc_deriv = cell2mat(rate_inc_deriv_cell);
%     thresh_deriv = floor(nanmedian(rate_inc_deriv));

    rate_inc_half_height_cell = cellfun(@(x,y) x<0.5*y,num2cell(single_co_smooth_up,2),num2cell(max_rate),'UniformOutput',false);
    rate_inc_half_height_cell = cellfun(@(x,y) find(x & y,1,'last'),rate_inc_half_height_cell,before_peak,'UniformOutput',false);
    problem_with_half_height{jj} = cellfun(@isempty,rate_inc_half_height_cell);
    rate_inc_half_height_cell(problem_with_half_height{jj}) = num2cell(NaN);
    rate_inc_half_height{jj} = cell2mat(rate_inc_half_height_cell);
    thresh_half_height{jj} = floor(nanmedian(rate_inc_half_height{jj}(~high_bl{jj})));
    
    mean_solo_cr{jj} = mean(solo_cr);
    ci_solo_cr = [prctile(solo_cr,solo_prcnt(1)),prctile(solo_cr,solo_prcnt(2))];
    err_solo{jj} = [ci_solo_cr(1)-mean_solo_cr{jj};mean_solo_cr{jj}-ci_solo_cr(2)];
    
    rate_inc_solo_cell = cellfun(@(x) x<ci_solo_cr(1),num2cell(single_co_smooth_up,2),'UniformOutput',false);
    rate_inc_solo_cell = cellfun(@(x,y) find(x & y,1,'last'),rate_inc_solo_cell,before_peak,'UniformOutput',false);
    problem_with_solo{jj} = cellfun(@isempty,rate_inc_solo_cell) | max_rate < ci_solo_cr(1);
    rate_inc_solo_cell(problem_with_solo{jj}) = num2cell(NaN);
    rate_inc_solo{jj} = cell2mat(rate_inc_solo_cell);
    thresh_solo{jj} = floor(nanmedian(rate_inc_solo{jj}(~high_bl{jj})));
    
    calling_rate(ii_dir).half_height_early{jj} = rate_inc_half_height{jj} < thresh_half_height{jj} & ~problem_with_half_height{jj};    
    calling_rate(ii_dir).half_height_late{jj} = rate_inc_half_height{jj} >= thresh_half_height{jj} & ~problem_with_half_height{jj};
    calling_rate(ii_dir).solo_early{jj} = rate_inc_solo{jj} < thresh_solo{jj} & ~problem_with_solo{jj};    
    calling_rate(ii_dir).solo_late{jj} = rate_inc_solo{jj} >= thresh_solo{jj} & ~problem_with_solo{jj};
    
    end
    
    %% single co figure
    click_dis = clicks_struct(1).co(ii_dir).clicks.dis_m;
    single_clicks = clicks_struct(1).co(ii_dir).clicks.not_DC;
    click_dis_single = click_dis .* single_clicks;
    %             click_dis(high_bl,:) = [];
    figure
    C = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
    for ico=1:length(max_rate)
        if high_bl{1}(ico)
            continue
        end
        clf
        axes('Position',[0.1 0.2 0.8 0.7])
        hold on
        xlabel('inter-bat distance (m)')
        %                 yyaxis left
        for jj=1:2
            p{jj}=shadedErrorBar([x_bins(1) x_bins(end)],repmat(mean_solo_cr{jj},1,2),repmat(err_solo{jj},1,2));
            set(p{jj}.edge,'visible','off')
            p{jj}.patch.FaceColor = C{jj};
            plot(x_bins,single_co_smooth{jj}(ico,:),'LineWidth',2,'color',C{jj})
            
            %                 yyaxis right
            %                 plot(x_bins(2:end),single_co_cr_diff(ico,:),'LineWidth',2,'color',[0.8500 0.3250 0.0980])
            %                 ylim([-7 7])
            %                 if ~problem_with_deriv(ico)
            %                     xline(x_bins(rate_inc_deriv(ico)),'--','color',[0.8500 0.3250 0.0980]);
            %                 end
            if ~problem_with_solo{jj}(ico)
                xline(x_bins_up(rate_inc_solo{jj}(ico)),'--','color',[0.5 0.5 0.5]);
            end
            if ~problem_with_half_height{jj}(ico)
                xline(x_bins_up(rate_inc_half_height{jj}(ico)),'--','color',C{jj});
            end
            if problem_with_solo{jj}(ico) && problem_with_half_height{jj}(ico)
                text(2,5,'problem with both threshold methods')
            end
        end
        xlim([x_bins(1) x_bins(end)])
        ylim([0 35])
        ylabel('Calling Rate (Hz)')
        %                 ylabel('Calling Rate diff (Hz)')
        title({name_prefix,['dir_' num2str(ii_dir) '_co_' num2str(ico)]},'Interpreter','none')
        legend(['Solo calling rate (Hz) - mean+' num2str(solo_prcnt(2)) '-' num2str(solo_prcnt(1)) '%'])
        
        axes('Position',[0.1 0.05 0.8 0.05])
        hold on
        plot([click_dis(ico,:);click_dis(ico,:)],repmat([1;-1],1,length(click_dis(ico,:))),'b');
        plot([click_dis_single(ico,:);click_dis_single(ico,:)],repmat([1;-1],1,length(click_dis_single(ico,:))),'r');
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
            thresh_used = thresh_half_height{1};
            co_early = calling_rate(ii_dir).half_height_early{1};
            co_late = calling_rate(ii_dir).half_height_late{1};
            rate_inc = rate_inc_half_height{1};
            criterion_used = 'half-height criterion';
            
            %                     thresh_used = thresh_deriv;
            %                     co_early = rate_inc_deriv < thresh_used;
            %                     rate_inc = rate_inc_deriv;
            %                     prob = problem_with_deriv;
            %                     criterion_used = ['diff criterion; diff thresh = ' num2str(diff_th)];
        else
            thresh_used = thresh_solo{1};
            co_early = calling_rate(ii_dir).solo_early{1};
            co_late = calling_rate(ii_dir).solo_late{1};
            rate_inc = rate_inc_solo{1};
            criterion_used = 'solo criterion';
        end
        hold on
        plot(x_bins,single_co_smooth{1}(co_early,:),'color',[1 0 0 .5],'LineWidth',0.7)
        plot(x_bins,single_co_smooth{1}(co_late,:),'color',[0 0 1 .5],'LineWidth',0.7)
        plot(x_bins,nanmean(single_co_smooth{1}(co_early,:)),'color',[1 0 0 1],'LineWidth',2.5)
        plot(x_bins,nanmean(single_co_smooth{1}(co_late,:)),'color',[0 0 1 1],'LineWidth',2.5)
        plot(x_bins_up(rate_inc(co_early)),cellfun(@(x,y) x(y),num2cell(single_co_smooth_up{1}(co_early,:),2),num2cell(rate_inc(co_early))),'ro')
        plot(x_bins_up(rate_inc(co_late)),cellfun(@(x,y) x(y),num2cell(single_co_smooth_up{1}(co_late,:),2),num2cell(rate_inc(co_late))),'bo')
        xline(x_bins_up(thresh_used),'--');
        xlabel('inter-bat distance (m)')
        ylabel('Calling Rate (Hz)')
        title({name_prefix,'CO events divided by calling-rate increase latency',['dir_' num2str(ii_dir) '; ' criterion_used]},'Interpreter','none')
        axes('Position',[.7 .7 .2 .2])
        hold on
        cla
%         [~,edges] = histcounts(x_bins_up(rate_inc(~isnan(rate_inc))));
%         histogram(x_bins_up(rate_inc(co_early)),edges,'FaceColor',[1 0 0])
%         histogram(x_bins_up(rate_inc(co_late)),edges,'FaceColor',[0 0 1])
        histogram(x_bins_up(rate_inc(co_early)),'BinWidth',2,'FaceColor',[1 0 0])
        histogram(x_bins_up(rate_inc(co_late)),'BinWidth',2,'FaceColor',[0 0 1])
        xlabel('inter-bat distance (m)')
        saveas(gcf,fullfile(click_rate_folder,['cr_hh_' num2str(i) '_' name_prefix,'_dir_' num2str(ii_dir) '.png']))
    end
    close
    
    %% click amplitude figure
    x_amp_bins = clicks_struct(1).co(ii_dir).calling_rate.intensity_dis_m_bins;
    
    solo_amp = clicks_struct(1).solo(ii_dir).click_intensity;
    median_solo_amp = nanmedian(solo_amp);
    err_solo = [prctile(solo_amp,solo_prcnt(1))-median_solo_amp;median_solo_amp-prctile(solo_amp,solo_prcnt(2))];
    
    co_amp = clicks_struct(1).co(ii_dir).calling_rate.single_co_intensity_dis_m;
    co_amp_median = nanmedian(co_amp,1);
    err_co = [prctile(co_amp,solo_prcnt(1),1)-co_amp_median;co_amp_median-prctile(co_amp,solo_prcnt(2),1)];
%     co_amp_mean = clicks_struct(1).co(ii_dir).calling_rate.intensity_dis_m;
%     co_amp_sm=smoothdata(co_amp,2,'movmean',smooth_wind);
%     co_amp_sm(isnan(co_amp)) = NaN;
    co_amp_median=smoothdata(co_amp_median,'movmean',smooth_wind);
    figure('WindowState','maximized')
    hold on
    p_solo=shadedErrorBar([x_bins(1) x_bins(end)],repmat(median_solo_amp,1,2),repmat(err_solo,1,2));
    set(get(get(p_solo.mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     p_co=plot(x_bins,co_amp_sm,'color',[0.3 0 0.3 .5],'LineWidth',0.7,'HandleVisibility','off');
%     plot(x_bins,co_amp_mean,'color',[0.3 0 0.3 1],'LineWidth',2.5)
    p_solo.mainLine.LineWidth = 2.5;
    set(p_solo.edge,'visible','off')
    p_co=shadedErrorBar(x_amp_bins,co_amp_median,err_co);
    p_co.patch.FaceColor = [0.3 0 0.3];
    set(p_co.mainLine,'LineWidth',2.5,'Color',[0.3 0 0.3])
    set(p_co.edge,'visible','off')
    xlabel('inter-bat distance (m)')
    ylabel('SNR')
    title({name_prefix,'Click Amplitude',['dir_' num2str(ii_dir)]},'Interpreter','none')
    L=legend([p_solo.mainLine,p_co.mainLine],{'Solo','CO'});
    title(L,['median + ' num2str(solo_prcnt(2)) '-' num2str(solo_prcnt(1)) '%'])
    saveas(gcf,fullfile(click_rate_folder,['amp_' name_prefix,'_dir_' num2str(ii_dir) '.png']))
    close
    
end