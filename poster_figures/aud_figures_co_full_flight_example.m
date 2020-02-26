clear
load('D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_audio_struct\bat_2336_20190902_CO_clicks_struct.mat')
%CO_struct=co_struct_new;
%%

figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [24 21]; % ~A4 size
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');
fsize=12;

self_color=[51 153 255]./255;
other_bat_color=[255 153 51]./255;
behav_color=[96 96 96;255 51 255]./255;

%%
for CO_i=1:length(CO_struct)
%behave_time_\\\limits=[1000    1.7541]*1.0e+09;
ts=CO_struct(CO_i).ts;
CO_point=CO_struct(CO_i).CO_point;
[~,ind_co_point]=min(abs(ts-CO_point));
lim_ind=[ind_co_point-6*10^5 ind_co_point+6*10^5];
ts=CO_struct(CO_i).ts(lim_ind(1):lim_ind(2));
behave_time_limits=CO_struct(CO_i).ts(lim_ind);
behave_time_ticks=linspace(behave_time_limits(1),behave_time_limits(2),13);
time_lim_str=round((behave_time_limits-mean(behave_time_limits))*(10^-6));
time_ticks_str=-6:6;
%% legend
y_legend=0.9;
axes('position',[0.75 y_legend 0.05 0.1]); 
plot([0 1], [0 0],'color',self_color,'LineWidth',3)
xlim([0 1])
text(0.5,0.3, 'Bat 1','HorizontalAlignment','center','FontWeight','bold','fontsize',fsize)
box off
axis off

axes('position',[0.81 y_legend 0.05 0.1]);
plot([0 1], [0 0],'color',other_bat_color,'LineWidth',3)
xlim([0 1])
text(0.5,0.3, 'Bat 2','HorizontalAlignment','center','FontWeight','bold','fontsize',fsize)
box off
axis off
%% pos
axes('position',[0.1 0.84 0.8 0.1])
plot(CO_struct(CO_i).bsp_ts,CO_struct(CO_i).pos_self,'color',self_color,'linewidth',3)
hold on;
plot(CO_struct(CO_i).bsp_ts,CO_struct(CO_i).pos_other,'color',other_bat_color,'linewidth',3)
plot([CO_point CO_point], [0 130],'LineStyle','--','color',behav_color(2,:),'linewidth', 2)
xlim(behave_time_limits)
ylim([10 130])
set(gca,'ytick',0:40:130)
set(gca,'xtick',[])
ylabel(sprintf('Position\nin tunnle (m)'))
box off;


%% aud large
axes('position',[0.1 0.62 0.8 0.2])
x=CO_struct(CO_i).filt_other(lim_ind(1):lim_ind(2));
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(ts ,xx,'color',other_bat_color)
hold on;
x=CO_struct(CO_i).filt_self(lim_ind(1):lim_ind(2));
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(ts,xx,'color',self_color)
plot([CO_point CO_point], [-1 1],'LineStyle','--','color',behav_color(2,:),'linewidth', 2)
set(gca,'xtick',behave_time_ticks,'xticklabel',time_ticks_str)
xlim(behave_time_limits)
ylim([-1 1])
ylabel('Audio signal (a.u)')
xlabel('Time (sec)')
box off;

%% aud close
distance=CO_struct(CO_i).distance ; 
lim_dis=[-30 10];

[~, lim_ind(1)]=min(abs(distance-lim_dis(1)));
[~, lim_ind(2)]=min(abs(distance-lim_dis(2)));

find_relevant_time=[CO_struct(CO_i).ts(lim_ind(1)) CO_struct(CO_i).ts(lim_ind(2))];
arrow_point=[CO_struct(CO_i).ts(lim_ind(1)),CO_struct(CO_i).ts(lim_ind(1));-1 -1];
[Xf,Yf]=ds2nfu(arrow_point(1,:),arrow_point(2,:));
Xf(2)=0.1;
Yf(2)=0.52;
annotation('line',Xf,Yf,'LineStyle','--','Color','k')

arrow_point=[CO_struct(CO_i).ts(lim_ind(2)),CO_struct(CO_i).ts(lim_ind(2));-1 -1];
[Xf,Yf]=ds2nfu(arrow_point(1,:),arrow_point(2,:));
Xf(2)=0.9;
Yf(2)=0.52;
annotation('line',Xf,Yf,'LineStyle','--','Color','k')

dis_to_plot=distance(lim_ind(1):lim_ind(2));
behave_dis_ticks=-30:10:10;

axes('position',[0.1 0.32 0.8 0.2])
x=CO_struct(CO_i).filt_other(lim_ind(1):lim_ind(2));
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(dis_to_plot ,xx,'color',other_bat_color)
hold on;
x=CO_struct(CO_i).filt_self(lim_ind(1):lim_ind(2));
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(dis_to_plot,xx,'color',self_color)
plot([0 0], [-1 1],'LineStyle','--','color',behav_color(2,:),'linewidth', 2)
set(gca,'xtick',behave_dis_ticks)
xlim(lim_dis)
ylim([-1 1])
ylabel('Audio signal (a.u)')
xlabel('Inter-bat distance (m)')
box off;

%% aud control
distance=CO_struct(CO_i).distance ; 
lim_dis=[-4 2];
behave_dis_ticks=-4:2:2;

[~, lim_ind(1)]=min(abs(distance-lim_dis(1)));
[~, lim_ind(2)]=min(abs(distance-lim_dis(2)));
dis_to_plot=distance(lim_ind(1):lim_ind(2));


arrow_point=[CO_struct(CO_i).distance(lim_ind(1)),CO_struct(CO_i).distance(lim_ind(1));-1 -1];
[Xf,Yf]=ds2nfu(arrow_point(1,:),arrow_point(2,:));
Xf(2)=0.1;
Yf(2)=0.23;
annotation('line',Xf,Yf,'LineStyle','--','Color','k')

arrow_point=[CO_struct(CO_i).distance(lim_ind(2)),CO_struct(CO_i).distance(lim_ind(2));-1 -1];
[Xf,Yf]=ds2nfu(arrow_point(1,:),arrow_point(2,:));
Xf(2)=0.9;
Yf(2)=0.23;
annotation('line',Xf,Yf,'LineStyle','--','Color','k')



bat_y=[0.07 0.16];
%plot
axes('position',[0.1 bat_y(1) 0.8 0.16])
clusters_other=CO_struct(CO_i).clusters_self;
clusters_self=CO_struct(CO_i).clusters_other;
x=CO_struct(CO_i).filt_other(lim_ind(1):lim_ind(2));
self_color_i=other_bat_color;
other_bat_color_i=self_color;
%peak of self:
find_relevant_ind=find([clusters_self.peak_abs_pos]>lim_ind(1) & [clusters_self.peak_abs_pos]<lim_ind(2));
relevat_peak_ind_self=[clusters_self.peak_abs_pos];
relevat_peak_ind_self=relevat_peak_ind_self(find_relevant_ind);
% top plot - self +peak of self
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(dis_to_plot,xx+1,'color',self_color_i);
hold on;
plot(distance(relevat_peak_ind_self),ones(1,length(relevat_peak_ind_self))+1,'.','color',self_color_i)

%bottom
x=CO_struct(CO_i).filt_self(lim_ind(1):lim_ind(2));
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(dis_to_plot,xx-1,'color',other_bat_color_i);

% peak of other:
find_relevant_ind=find([clusters_other.peak_abs_pos]>lim_ind(1) & [clusters_other.peak_abs_pos]<lim_ind(2));
relevat_peak_ind_self=[clusters_other.peak_abs_pos];
relevat_peak_ind_self=relevat_peak_ind_self(find_relevant_ind);

%find_relevant_ind=find([clusters_other.abs_pos_at_2nd_bat]>lim_ind(1) & [clusters_other.abs_pos_at_2nd_bat]<lim_ind(2));
relevat_peak_ind_other=[clusters_other.abs_pos_at_2nd_bat];
relevat_peak_ind_other=relevat_peak_ind_other(find_relevant_ind);

for ii=1:length(relevat_peak_ind_self)
   plot([distance(relevat_peak_ind_self(ii)) distance(relevat_peak_ind_other(ii))],[-0.2 0.2],'k.-') 
end
plot([0 0], [-2.2 2.2],'LineStyle','--','color',behav_color(2,:),'linewidth', 2)
xlabel('Inter-bat distance (m)')
xlim(lim_dis)
ylim([-2.2 2.2])
set(gca,'ytick',[])
set(gca,'xtick',behave_dis_ticks)
box off;
% 
% for bat_i=1:2
%     if bat_i==1
%         clusters_self=CO_struct(CO_i).clusters_self;
%         clusters_other=CO_struct(CO_i).clusters_other;
%         x=CO_struct(CO_i).filt_self(lim_ind(1):lim_ind(2));
%         self_color_i=self_color;
%         other_bat_color_i=other_bat_color;
%     else
%         clusters_other=CO_struct(CO_i).clusters_self;
%         clusters_self=CO_struct(CO_i).clusters_other;
%         x=CO_struct(CO_i).filt_other(lim_ind(1):lim_ind(2));
%         self_color_i=other_bat_color;
%         other_bat_color_i=self_color;
%     end
%     %peak of self:
% find_relevant_ind=find([clusters_self.peak_abs_pos]>lim_ind(1) & [clusters_self.peak_abs_pos]<lim_ind(2));
% relevat_peak_ind_self=[clusters_self.peak_abs_pos];
% relevat_peak_ind_self=relevat_peak_ind_self(find_relevant_ind);
% 
% % peak of other:
% find_relevant_ind=find([clusters_other.abs_pos_at_2nd_bat]>lim_ind(1) & [clusters_other.abs_pos_at_2nd_bat]<lim_ind(2));
% relevat_peak_ind_other=[clusters_other.abs_pos_at_2nd_bat];
% relevat_peak_ind_other=relevat_peak_ind_other(find_relevant_ind);
% 
% 
% %plot
% axes('position',[0.1 bat_y(bat_i) 0.8 0.08])
% hold on;
% xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
% xx=xx-nanmean(xx);
% plot(dis_to_plot,xx,'color',self_color_i)
% plot(distance(relevat_peak_ind_self),ones(1,length(relevat_peak_ind_self)),'.','color',self_color_i)
% for ii=1:length(relevat_peak_ind_other)
% plot([distance(relevat_peak_ind_other(ii)) distance(relevat_peak_ind_other(ii))],[-1 1],'--','color',other_bat_color_i)
% end
% plot([0 0], [-1 1],'LineStyle','--','color',behav_color(2,:),'linewidth', 2)
% xlim(lim_dis)
% ylim([-1.1 1.1])
% ylabel('Audio signal (a.u)')
% if bat_i==1
%     set(gca,'xtick',behave_dis_ticks)
% xlabel('Inter-bat distance (m)')
% else
%         set(gca,'xtick',[])
% 
% end
% box off;
% 
% end


%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,['aud_example_co',num2str(CO_i)]);
print(gcf, fig_name, '-dpdf',  '-painters');

clf   
end