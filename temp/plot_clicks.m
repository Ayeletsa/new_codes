function plot_clicks(CO_struct,clusters_self,clusters_other,points2plot)
hold on
plot(CO_struct.ts,CO_struct.filt_self,'r')
plot(CO_struct.ts,CO_struct.filt_other,'g')
switch points2plot
    case 'max'
        plot(CO_struct.ts([clusters_self.peak_abs_pos]),[clusters_self.MaxIntensity],'*k')
        plot(CO_struct.ts([clusters_other.peak_abs_pos]),[clusters_other.MaxIntensity],'*k')
        C_s=num2cell(1:numel(clusters_self));
        dx = 0.2; dy = 0.2;
        text(CO_struct.ts([clusters_self.peak_abs_pos])+dx,[clusters_self.MaxIntensity]+dy,C_s);
        C_o=num2cell(1:numel(clusters_other));
        text(CO_struct.ts([clusters_other.peak_abs_pos])+dx,[clusters_other.MaxIntensity]+dy,C_o);
    case 'all'
        plot(CO_struct.ts(cat(1,clusters_self.PixelIdxList)),cat(2,clusters_self.PixelValues),'*k')
        plot(CO_struct.ts(cat(1,clusters_other.PixelIdxList)),cat(2,clusters_other.PixelValues),'*k')
    case 'align'
        plot(CO_struct.ts([clusters_self.peak_abs_pos]),[clusters_self.MaxIntensity],'*m')
        plot(CO_struct.ts([clusters_other.peak_abs_pos]),[clusters_other.MaxIntensity],'*b')
%         plot(CO_struct.ts([clusters_self.abs_pos_at_2nd_bat]),[clusters_self.MaxIntensity],'Xb')
%         plot(CO_struct.ts([clusters_other.abs_pos_at_2nd_bat]),[clusters_other.MaxIntensity],'Xm')
        plot(CO_struct.ts([clusters_self.abs_pos_at_2nd_bat;clusters_self.abs_pos_at_2nd_bat]),[zeros(1,length(clusters_self));clusters_self.MaxIntensity],'--b')
        plot(CO_struct.ts([clusters_other.abs_pos_at_2nd_bat;clusters_other.abs_pos_at_2nd_bat]),[zeros(1,length(clusters_other));clusters_other.MaxIntensity],'--m')
    otherwise
end
title(['CO #' num2str(CO_struct.num)])
end