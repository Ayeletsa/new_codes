clear

file_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_audio_struct\CO_audio_population.mat';
load(file_name)

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
colors=[self_color;other_bat_color];
%%
bin_centers=H(1).bin_centers  ;
for ii=2:length(H)
    aud(ii-1,:,2)=H(ii).bat31_smooth;
    aud(ii-1,:,1)=H(ii).bat36_smooth;
    
end

%%
y_pos=[0.5 0.2];
for bat_i=1:2
    aud_bat=aud(:,:,bat_i);
    mean_aud=mean(aud_bat);
    
    axes('position',[0.1 y_pos(bat_i) 0.2 0.2]);
    plot(bin_centers,aud_bat,'color',colors(bat_i,:),'linewidth',2)
    hold on;
    plot(bin_centers,mean_aud,'color','k','linewidth',1)
    xlim([-40 40])
    xlabel('Inter-bat distance (m)')
    ylabel('Click rate (Hz)')
    
end

%% save figure

fig = gcf;
%fig.InvertHardcopy = 'off';
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end
fig_name=fullfile(dir_name,['aud_population']);
print(gcf, fig_name, '-dpdf',  '-painters');

clf  
