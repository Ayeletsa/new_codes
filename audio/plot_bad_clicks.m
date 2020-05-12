if ~exist('fig_prefix','var')
    fig_prefix = clicks_struct(1).name;
end

% plot examples of bad clicks
plot_margin = 80;
max_per_fig = 25;
aspect_ratio = 1/1;

for ibat=1:2
%     clicks_to_use = clicks_struct(ibat).wrong_rise_time;
%     fig_suffix = 'rise_time';
%     plot_many_clicks(clicks,ibat,clicks_struct,clicks_to_use,plot_margin,max_per_fig,aspect_ratio,fullfile(fig_folder,[fig_prefix '_' clicks(ibat).bat '_' fig_suffix]))
%     
%     clicks_to_use = clicks_struct(ibat).wrong_size;
%     fig_suffix = 'length';
%     plot_many_clicks(clicks,ibat,clicks_struct,clicks_to_use,plot_margin,max_per_fig,aspect_ratio,fullfile(fig_folder,[fig_prefix '_' clicks(ibat).bat '_' fig_suffix]))
%     
%     clicks_to_use = clicks_struct(ibat).wrong_spectrum;
%     fig_suffix = 'spectrum';
%     plot_many_clicks(clicks,ibat,clicks_struct,clicks_to_use,plot_margin,max_per_fig,aspect_ratio,fullfile(fig_folder,[fig_prefix '_' clicks(ibat).bat '_' fig_suffix]))

end

function plot_many_clicks(clicks,ibat,clicks_struct,clicks_to_use,plot_margin,max_per_fig,aspect_ratio,fig_name,plot_spectrum)
grc = @(x,r) [ceil(sqrt(r*x)),ceil(sqrt(r*x)/r)];
rc = grc(max_per_fig,aspect_ratio);
max_per_fig = prod(rc);
ncl = length(clicks_to_use);
q = floor(ncl/max_per_fig);
r = rem(ncl,max_per_fig);
v = [repmat(max_per_fig,1,q),r];
ve = cumsum(v);
vs = ve - v + 1;

figure('WindowState','maximized')
for ifig=1:length(v)
    clf
    cl_vec = vs(ifig):ve(ifig);
    for icl=1:length(cl_vec)
        cl = cl_vec(icl);
        subplot(rc(1),rc(2),icl);
        hold on
        cl_p = clicks_to_use(cl).peak_abs_idx;
        cl_s = clicks_to_use(cl).start;
        cl_e = clicks_to_use(cl).end;
        plot_wind = cl_s - plot_margin : cl_e + plot_margin;
        t0 = plot_wind(1);

        plot((plot_wind-t0)*1e-2,clicks(ibat).filt(plot_wind)/clicks_struct(ibat).aud_BL)
        plot((cl_s-t0)*1e-2,clicks(ibat).filt(cl_s)/clicks_struct(ibat).aud_BL,'ro')
        plot((cl_p-t0)*1e-2,clicks(ibat).filt(cl_p)/clicks_struct(ibat).aud_BL,'ro')
        plot((cl_e-t0)*1e-2,clicks(ibat).filt(cl_e)/clicks_struct(ibat).aud_BL,'ro')
        
%         ylim([-200 200])
        title(num2str(cl))
%         gca.XAxis.Visible = 'off';
        if icl==length(cl_vec)
            xlabel('ms')
            ylabel('snr')
        end
        
    end
    [~,t,~] = fileparts(fig_name);
    sgtitle(t,'interpreter','none')
    saveas(gcf,[fig_name '_' num2str(ifig) '.png'])
end
close
end
