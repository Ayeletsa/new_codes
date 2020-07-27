
function data=PV_get_cell_data(cell_co_solo_initial_analysis)
for dir_i=1:2
%load flight solo data:
%bsp:
data(dir_i).solo.bsp.flight_ts=cell_co_solo_initial_analysis.solo(dir_i).bsp.ts_usec;
data(dir_i).solo.bsp.flight_x_pos=cell_co_solo_initial_analysis.solo(dir_i).bsp.x_pos;
data(dir_i).solo.bsp.flight_y_pos=cell_co_solo_initial_analysis.solo(dir_i).bsp.y_pos;
%spikes:
data(dir_i).solo.spikes.flight_x_pos=cell_co_solo_initial_analysis.solo(dir_i).spikes.x_pos;
data(dir_i).solo.spikes.flight_y_pos=cell_co_solo_initial_analysis.solo(dir_i).spikes.y_pos;
data(dir_i).solo.spikes.flight_ts=cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec;

%load co data:
%bsp:
data(dir_i).co.bsp.flight_ts=cell_co_solo_initial_analysis.co(dir_i).bsp_full.ts_usec;
data(dir_i).co.bsp.flight_x_pos=cell_co_solo_initial_analysis.co(dir_i).bsp_full.x_pos;
data(dir_i).co.bsp.flight_y_pos=cell_co_solo_initial_analysis.co(dir_i).bsp_full.y_pos;
data(dir_i).co.bsp.flight_dis=cell_co_solo_initial_analysis.co(dir_i).bsp_full.dis_m;

%spikes:
data(dir_i).co.spikes.flight_ts=cell_co_solo_initial_analysis.co(dir_i).spikes_full.ts_usec;
data(dir_i).co.spikes.flight_x_pos=cell_co_solo_initial_analysis.co(dir_i).spikes_full.x_pos;
data(dir_i).co.spikes.flight_y_pos=cell_co_solo_initial_analysis.co(dir_i).spikes_full.y_pos;
data(dir_i).co.spikes.flight_dis=cell_co_solo_initial_analysis.co(dir_i).spikes_full.dis_m;

% velocity short distance prop to medien
data(dir_i).co.bsp.vel_prop_to_median=cell_co_solo_initial_analysis.co(dir_i).info.vel.vel_prop_to_median ; 
end