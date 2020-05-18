function vel=compute_velocity(x,y,ts,factor_to_sec)

vel=((sqrt(diff(x).^2+diff(y).^2)./diff(ts)))*factor_to_sec;
