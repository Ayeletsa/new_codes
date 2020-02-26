clear
figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [30 21]; 
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');
fsize=14;
%behavior analysis fig for sfn poster
example_day_modes='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\behavioral_modes_bat_2389_day_20180806.mat';
example_day_general='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\general_behave_analysis_bat2389_day20180806.mat';
load(example_day_modes)
load(example_day_general)
%relevant times for plot:
FE_ind=bsp_proc_data(tag_i).flight_ind;
ts=bsp_proc_data(tag_i).ts;
number_of_plots=9;
end_time=3.0070*1e10;
total_time=end_time-ts(FE_ind(1));
time_bin=total_time/number_of_plots;
time_vec_for_plot=ts(FE_ind(1)):time_bin:ts(FE_ind(end));
%% relevant behavior


relevant_behavior_str={'Solo flights','Cross-overs'};
[~,behav_modes_plot(1).start,behav_modes_plot(1).end]=find_length_of_consecutive_ind(behavioral_modes.solo_ind,length(pos_self_x));
[~,behav_modes_plot(2).start,behav_modes_plot(2).end]=find_length_of_consecutive_ind(behavioral_modes.CO_ind,length(pos_self_x));

%% figure param:
panel_size=[0.8 0.09];
x_position=0.05;
behav_color=[96 96 96;255 51 255]./255;
self_color=[51 153 255]./255;
other_bat_color=[255 153 51]./255;

vertical_dist=0.095;
y_position(1)=0.85;
y_position(2)=y_position(1)-vertical_dist;
y_position(3)=y_position(2)-vertical_dist;
y_position(4)=y_position(3)-vertical_dist;
y_position(5)=y_position(4)-vertical_dist;
y_position(6)=y_position(5)-vertical_dist;
y_position(7)=y_position(6)-vertical_dist;
y_position(8)=y_position(7)-vertical_dist;
y_position(9)=y_position(8)-vertical_dist;

%% PLOT
%figure('units','normalized','outerposition',[0 0 1 1])

% if tag_i==2

% title:
% axes('position', [0.4, y_position(1), 0.1,  0.1]);
% %title(sprintf('Bat %d day %d - Total recording time %.1f min, Total flight time %.1f min',bat,day,(total_time/us_factor)/60,(length(FE_ind)/frame_per_second)/60), 'FontSize', 11, 'FontWeight', 'Bold');
% box off;
% axis off;
% else
%     color_self=[1 0 0];
%     color_other=[0 1 0];
% end
% All behavior
%--------------------------------------------------

for time_i=1:length(time_vec_for_plot)-1
    ax(time_i)=axes('position',[x_position y_position(time_i) panel_size]);
    relevant_behav_ind=find(ts>=time_vec_for_plot(time_i) & ts<time_vec_for_plot(time_i+1));
    by_pass=behavioral_modes.bypass_ind(find(ts(behavioral_modes.bypass_ind)>=time_vec_for_plot(time_i) & ts(behavioral_modes.bypass_ind)<time_vec_for_plot(time_i+1)));
    plot(ts(relevant_behav_ind),pos_other_x(relevant_behav_ind),'.','color',other_bat_color)
    hold on;
    plot(ts(relevant_behav_ind),pos_self_x(relevant_behav_ind),'.','color',self_color);
    hold on;
    if sum(isnan(pos_self_x(relevant_behav_ind)))==length(pos_self_x(relevant_behav_ind))
        set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))])
    else
        set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))],'ylim' ,[min(pos_self_x(relevant_behav_ind)) max(pos_self_x(relevant_behav_ind))])
    end
    set(gca,'XTick',[min(ts) max(ts)],'XTickLabel',round(([min(ts) max(ts)])/1e3/60))
    ylimits=[min([pos_self_x; pos_other_x]) max([pos_self_x; pos_other_x])];
    ylim(ylimits)
   set(gca,'YTick',[])

    if ~isempty(by_pass)
        plot(ts(by_pass),pos_self_x(by_pass),'*k');
    end
    
    behav_mod_i=1;
        for ii=1:length(behav_modes_plot(behav_mod_i).start)
            start_point=[];
            end_point=[];
            if ts(behav_modes_plot(behav_mod_i).start(ii))>=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))<=time_vec_for_plot(time_i+1)
                start_point=ts(behav_modes_plot(behav_mod_i).start(ii));
                end_point=ts(behav_modes_plot(behav_mod_i).end(ii));
            elseif ts(behav_modes_plot(behav_mod_i).start(ii))>=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))>=time_vec_for_plot(time_i+1)
                start_point=ts(behav_modes_plot(behav_mod_i).start(ii));
                end_point=time_vec_for_plot(time_i+1);
            elseif ts(behav_modes_plot(behav_mod_i).start(ii))<=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))<=time_vec_for_plot(time_i+1)
                start_point=time_vec_for_plot(time_i);
                end_point=ts(behav_modes_plot(behav_mod_i).end(ii));
                
                
            end
            if ~isempty(start_point)
                rec_x=[start_point  end_point  end_point start_point];
                rec_y=[min(pos_self_x)  min(pos_self_x) max(pos_self_x) max(pos_self_x)];
                p=patch(rec_x,rec_y,behav_color(behav_mod_i,:),'EdgeColor','none');
                set(p,'FaceAlpha',0.3)
            end
        % CO point
        plot(ts(behavioral_modes.CO_point), pos_self_x(behavioral_modes.CO_point),'o','markersize',8,'MarkerFaceColor',behav_color(2,:),'MarkerEdgeColor','k')
        
    end
end

%% create time scale:
axes('position',[x_position 0.07 panel_size(1) 0.02]);
time_bin_in_min=(time_bin/us_factor)/60;
plot([0 time_bin/time_bin_in_min], [0 0],'k','LineWidth',5)
xlim([0 time_bin])
text(time_bin/(2*time_bin_in_min),-1.5, '1 min','HorizontalAlignment','center','fontsize',fsize*0.8)
box off
axis off

%% create legend
y_legend=0.91;

axes('position',[0.57 y_legend 0.05 0.1]);
plot([0 1], [0 0],'color',self_color,'LineWidth',3)
xlim([0 1])
text(0.5,0.3, 'Bat 1','HorizontalAlignment','center','FontWeight','bold','fontsize',fsize)
box off
axis off

axes('position',[0.63 y_legend 0.05 0.1]);
plot([0 1], [0 0],'color',other_bat_color,'LineWidth',3)
xlim([0 1])
text(0.5,0.3, 'Bat 2','HorizontalAlignment','center','FontWeight','bold','fontsize',fsize)
box off
axis off

axes('position',[0.69 y_legend+0.05 0.05 0.03]);
rec_x=[0,1,1,0];
rec_y=[0,0,1,1];
p=patch(rec_x,rec_y,behav_color(1,:),'EdgeColor','none');
set(p,'FaceAlpha',0.3)
text(0.5,0.5, 'Solo','HorizontalAlignment','center','color', [0 0 0],'FontWeight','bold','fontsize',fsize)
box off
axis off

axes('position',[0.72 y_legend+0.045 0.07 0.03]);
plot(0.1,0.5,'o','markersize',10,'MarkerFaceColor',behav_color(2,:),'MarkerEdgeColor','k')
text(0.5,0.5, 'Cross-overs','color', [0 0 0],'FontWeight','bold','fontsize',fsize)
box off
axis off



%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,'behav_example');
print(gcf, fig_name, '-dpdf',  '-painters');

clf
