clear
load('D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_audio_struct\bat_2336_20190902_CO_clicks_struct_less_co.mat')
CO_struct=co_struct_new;
%%
figure('units','normalized','outerposition',[0 0 1 1])
 folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\co_aud';
for CO_i=3
ax1=subplot(2,1,1);
x=CO_struct(CO_i).filt_other;
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(CO_struct(CO_i).ts ,xx,'r')
hold on;
x=CO_struct(CO_i).filt_self;
xx=(2./(nanmax(x)-nanmin(x))).*(x-nanmin(x))-1;
xx=xx-nanmean(xx);
plot(CO_struct(CO_i).ts,xx,'g')
xlim([min(CO_struct(CO_i).ts) max(CO_struct(CO_i).ts)])

ax2=subplot(2,1,2);
plot(CO_struct(CO_i).bsp_ts,CO_struct(CO_i).pos_other,'r')
hold on;
plot(CO_struct(CO_i).bsp_ts,CO_struct(CO_i).pos_self,'g')
xlim([min(CO_struct(CO_i).ts) max(CO_struct(CO_i).ts)])

linkaxes([ax1,ax2],'x')

fig_name=fullfile(folder,['CO_',num2str(CO_i),'.jpg']);
saveas(gcf,fig_name)
clf
    
end