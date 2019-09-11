figure('units','normalized','outerposition',[0 0 1 1])

ax1=subplot(2,1,1);
plot(p.Aud.self_ts(end/8:end/4),p.Aud.self_signal(end/8:end/4),'r')
hold on;
plot(p.Aud.other_ts_us(end/8:end/4),p.Aud.other_signal(end/8:end/4),'g')
xlim([min(p.Aud.self_ts(end/8:end/4)) max(p.Aud.self_ts(end/8:end/4))])

ax2=subplot(2,1,2);
tag_self=find(p.bsp_tag_self==[p.bsp_data.tag_ID]);

plot(p.bsp_data(tag_self).ts_us_upsampled,p.bsp_data(tag_self).pos_upsampled(1,:),'r')
hold on;

tag_other=find(p.bsp_tag_other==[p.bsp_data.tag_ID]);

plot(p.bsp_data(tag_other).ts_us_upsampled,p.bsp_data(tag_other).pos_upsampled(1,:),'g')
xlim([min(p.Aud.self_ts(end/8:end/4)) max(p.Aud.self_ts(end/8:end/4))])

linkaxes([ax1,ax2],'x')
