function [behavioral_modes]=find_behavioral_modes(bsp_data,behav_params_file_name,tag_i,bat,day,behave_ts,dir_param_file_name)
%load parameters:
load(behav_params_file_name)
load(dir_param_file_name)
%% figure param:
panel_size=[0.92 0.095];
x_position=0.02;
direc_color=[0 0.447 0.741;0.8500    0.3250    0.0980];

vertical_dist=0.1;
y_position(1)=0.85;
y_position(2)=y_position(1)-vertical_dist;
y_position(3)=y_position(2)-vertical_dist;
y_position(4)=y_position(3)-vertical_dist;
y_position(5)=y_position(4)-vertical_dist;
y_position(6)=y_position(5)-vertical_dist;
y_position(7)=y_position(6)-vertical_dist;
y_position(8)=y_position(7)-vertical_dist;
y_position(9)=y_position(8)-vertical_dist;

number_of_plots=9;
us_factor=1e6;
%% define general behavioral variables:

FE_ind=bsp_data(tag_i).flight_ind;
ts=bsp_data(tag_i).ts;


%relevant times for plot:
total_time=ts(FE_ind(end))-ts(FE_ind(1));
time_bin=total_time/number_of_plots;
time_vec_for_plot=ts(FE_ind(1)):time_bin:ts(FE_ind(end));

%flights position and direction self:
%-----------------------------------------------------------------------------
only_this_bat_is_flying_ind_FE=setdiff(FE_ind,bsp_data(3-tag_i).flight_ind);
only_this_bat_is_flying_ind_FE=find(ismember(FE_ind,only_this_bat_is_flying_ind_FE));
%find the bat X position and velocity
pos_self_x=bsp_data(tag_i).pos(:,1);
pos_self_y=bsp_data(tag_i).pos(:,2);
pos_self_x_FE=pos_self_x(FE_ind);
velocity_self=(diff([0; pos_self_x])./diff([0;bsp_data(tag_i).ts])).*us_factor;
velocity_self(1)=NaN;
velocity_self(isnan(pos_self_x))=NaN;
velocity_self_FE=velocity_self(FE_ind);
% find the xy velocity
pos_self_xy=bsp_data(tag_i).pos(:,1:2);
velocity_self_xy=(sqrt(nansum(diff([0 0;pos_self_xy]).^2,2))./diff([0;bsp_data(tag_i).ts])).*us_factor;
velocity_self_xy(1)=NaN;
velocity_self_xy(isnan(pos_self_x))=NaN;

% find flight directions:
behavioral_modes.directional_ind{1}=FE_ind(find(velocity_self_FE>0));
behavioral_modes.directional_ind{2}=FE_ind(find(velocity_self_FE<0));

%flights position and direction other:
%-----------------------------------------------------------------------------
pos_other_x=bsp_data(3-tag_i).pos(:,1);
pos_other_y=bsp_data(3-tag_i).pos(:,2);

pos_other_x_FE=pos_other_x(FE_ind);
velocity_other=(diff([0; bsp_data(3-tag_i).pos(:,1)])./diff([0;bsp_data(3-tag_i).ts])).*us_factor;
velocity_other(1)=NaN;
velocity_other_FE=velocity_other(FE_ind);
velocity_other_FE(only_this_bat_is_flying_ind_FE)=0;
velocity_other(isnan(pos_other_x))=NaN;

%Distance between the bats:
%-----------------------------------------------------------------------------
distnace_other_from_self_FE=pos_other_x_FE-pos_self_x_FE;
distnace_other_from_self=pos_other_x-pos_self_x;

distnace_other_from_self_y=pos_other_y-pos_self_y;

behavioral_modes.pos_other=bsp_data(3-tag_i).pos;
behavioral_modes.pos_self=bsp_data(tag_i).pos;
behavioral_modes.ts=ts;
%% Solo
%Solo is defined when the bats is flying alone or when both are flying but they
%are far away from each other:
potential_solo_ind=FE_ind(unique([find(abs(distnace_other_from_self_FE)>dist_thresh_solo);only_this_bat_is_flying_ind_FE]));
[solo_length,solo_start,solo_end]=find_length_of_consecutive_ind(potential_solo_ind,length(pos_self_x));
long_solo_ind=find(solo_length>min_solo_length);
% take only long flights:
solo_ind=[];
for track_i=1:length(long_solo_ind)
    solo_ind=[solo_ind, solo_start(long_solo_ind(track_i)):solo_end(long_solo_ind(track_i))];
end

behavioral_modes.solo_ind=solo_ind;

behav=1;
behav_modes_plot(behav).name='Solo';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(solo_ind,length(pos_self_x));

%% Tracking
% tracking is defined when the two bats are flying and they are in the same direction (velocity
% with the same sign) and the distance between them is bellow threshold for
% more than t sec

%Conditions:
%1) Both bats are flying:
both_bats_are_flying=intersect(FE_ind,bsp_data(3-tag_i).flight_ind);
both_bats_are_flying_FE=find(ismember(FE_ind,both_bats_are_flying));
%2) short distance:
short_distance_ind_FE=find(abs(distnace_other_from_self_FE)<dist_thresh_tracking);
%3) same flight direction:
same_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==1);

potential_tracking_ind_FE=intersect(same_flight_direc_ind_FE,both_bats_are_flying_FE);
potential_tracking_ind_FE=intersect(potential_tracking_ind_FE,short_distance_ind_FE);

%4) divide to when the bat is tracking of being tracked:
%a. relative distance between the bats:
bat_behind=FE_ind([find(distnace_other_from_self_FE>0 & velocity_self_FE>0);find(distnace_other_from_self_FE<0 & velocity_self_FE<0)]);
bat_ahead=FE_ind([find(distnace_other_from_self_FE<0 & velocity_self_FE>0);find(distnace_other_from_self_FE>0 & velocity_self_FE<0)]);
%b. combine with previous conditions:
being_tracked_ind_also_short=intersect(FE_ind(potential_tracking_ind_FE),bat_ahead);
tracking_other_ind_also_short=intersect(FE_ind(potential_tracking_ind_FE),bat_behind);

%5) minimum time for tracking:
%initialize to find start and end of tracking:
%a. being tracked:
%-------------------------------
[being_tracked_length,being_tracked_start,being_tracked_end]=find_length_of_consecutive_ind(being_tracked_ind_also_short,length(pos_self_x));
long_being_tracked_ind=find(being_tracked_length>min_tracking_length);
being_tracked_ind=[];
for track_i=1:length(long_being_tracked_ind)
    being_tracked_ind=[being_tracked_ind, being_tracked_start(long_being_tracked_ind(track_i)):being_tracked_end(long_being_tracked_ind(track_i))];
end

%b. tracking
%-------------------------------
[tracking_other_length,tracking_other_start,tracking_other_end]=find_length_of_consecutive_ind(tracking_other_ind_also_short,length(pos_self_x));
long_tracking_other_ind=find(tracking_other_length>min_tracking_length);
tracking_other_ind=[];
for track_i=1:length(long_tracking_other_ind)
    tracking_other_ind=[tracking_other_ind, tracking_other_start(long_tracking_other_ind(track_i)):tracking_other_end(long_tracking_other_ind(track_i))];
end

% assign to structs:
behavioral_modes.being_tracked_ind=being_tracked_ind;
behavioral_modes.tracking_other_ind=tracking_other_ind;
%combine:
behavioral_modes.tracking_ind=sort([tracking_other_ind,being_tracked_ind]);
tracking_ind=behavioral_modes.tracking_ind;

% Find ind of bypass-
%-------------------------------------------------------------------------
% Defined by change from tracking to being tracked

%a. find where distance change signs:
distnace_other_from_self_shifted_plus=[0; distnace_other_from_self];
distnace_other_from_self_shifted_minus=[distnace_other_from_self; 0];
distance_change_sign=find(sign(distnace_other_from_self_shifted_minus.*distnace_other_from_self_shifted_plus)==-1)-1;
% find intersections with tracking:
behavioral_modes.bypass_ind=intersect(distance_change_sign,tracking_ind);

%find start_end for plots:
behav=2;
behav_modes_plot(behav).name='tracking';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(behavioral_modes.tracking_other_ind,length(pos_self_x));
behav=3;
behav_modes_plot(behav).name='being tracked';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(behavioral_modes.being_tracked_ind,length(pos_self_x));

%% Cross over:
%Cross over is defined when both bats are flying in differnt direction for t time before and t time after
%the distance change sign and te bats are in different directions:

%1) Both bats are flying (found above):
%2) bats with different flight directions:
different_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==-1);
[different_direc_length,different_direc_start,different_direc_end]=find_length_of_consecutive_ind(FE_ind(different_flight_direc_ind_FE),length(pos_self_x));

%3) intercet prev conditions:
relevant_flights_FE=intersect(both_bats_are_flying_FE,different_flight_direc_ind_FE);
%relevant_flights_FE=intersect(relevant_flights_FE,short_distance_ind_FE);

%4) distance change sign:
distnace_other_from_self_shifted_plusFE=[0; distnace_other_from_self_FE];
distnace_other_from_self_shifted_minusFE=[distnace_other_from_self_FE; 0];

distance_change_sign=FE_ind(find(sign(distnace_other_from_self_shifted_minusFE.*distnace_other_from_self_shifted_plusFE)==-1)-1);
CO_ind=[];
CO_point=[];

count=0;
for CO_i=1:length(distance_change_sign)
    
    relevant_event=find(distance_change_sign(CO_i)>different_direc_start & distance_change_sign(CO_i)<different_direc_end);
    if ~isempty(relevant_event)
        %find if before and after the CO they are in different directions(no UT)
        differnt_direc_before_CO=distance_change_sign(CO_i)-different_direc_start(relevant_event);
        differnt_direc_after_CO=different_direc_end(relevant_event)-distance_change_sign(CO_i);
        if differnt_direc_before_CO>min_time_before_CO & differnt_direc_after_CO>min_time_before_CO
            count=count+1;
            CO_point=[CO_point distance_change_sign(CO_i)];
            %find window:
            % find the start:
            long_dist_ind=find(abs(distnace_other_from_self)>=CO_window(1));
            before_CO=long_dist_ind(find([long_dist_ind-CO_point(end)]<0));
            start_ind=before_CO(end);%first ind before CO that the distance between the bats is lrager than ..
            if (CO_point(count)-start_ind)>max_wind_CO
                start_ind=CO_point(count)-max_wind_CO;
            end
            %find the end of the window:
            long_dist_ind=find(abs(distnace_other_from_self)>=CO_window(2));
            after_CO=long_dist_ind(find([long_dist_ind-CO_point(end)]>0));
            if ~isempty(after_CO)
                end_ind=after_CO(1); %first ind after CO that the distance between the bats is lrager than ..
                
                CO_ind=[CO_ind start_ind:end_ind];
                CO_event{CO_i}=start_ind:end_ind;
            end
        end
    end
end
%remove empty cells:
CO_event(cellfun('isempty',CO_event))=[];

behavioral_modes.CO_event=CO_event;
behavioral_modes.CO_ind=CO_ind;
behavioral_modes.CO_point=CO_point;
%% U-turns
% U-turns are defined like CO that one of the bat changed his direction X
% meters after CO
%1. find ind that bats change directions (vel changed sign):
% %self:
velocity_self_shifted_plus=[velocity_self; 0];
velocity_self_shifted_minus=[0 ;velocity_self];
change_direction_self_ind=find(sign(velocity_self_shifted_plus.*velocity_self_shifted_minus)==-1)-1;
% %other:
velocity_other_shifted_plus=[velocity_other; 0];
velocity_other_shifted_minus=[0 ;velocity_other];
change_direction_other_ind=find(sign(velocity_other_shifted_plus.*velocity_other_shifted_minus)==-1)-1;


UT_ind=[];
UT_point=[];
differnt_direc_before_CO=[];
CO_point_before_UT=[];
count=0;
ii=0;
iii=0;
for UT_i=1:length(distance_change_sign)
    
    relevant_event=find(distance_change_sign(UT_i)>different_direc_start & distance_change_sign(UT_i)<different_direc_end);
    if ~isempty(relevant_event)
        %find if before the CO they are in different directions(no UT before)
        differnt_direc_before_CO=distance_change_sign(UT_i)-different_direc_start(relevant_event);
        %         [closest_UT_ts, ind]=min(abs(change_direction_self_ind-distance_change_sign(CO_i));
        %         closest_UT_pos=abs(pos_self(ind)-distance_change_sign(CO_i));
        
        %1. test if they were in different direction before CO:
        if differnt_direc_before_CO>min_time_before_CO
            CO_point_before_UT(UT_i)=distance_change_sign(UT_i);
            count=count+1;
            
            %UT self:
            %1. find when in time there were UTs
            
            time_after_CO_from_self_CD=change_direction_self_ind-CO_point_before_UT(UT_i);
            time_after_CO_from_self_CD(time_after_CO_from_self_CD<0)=NaN; %direction changed before CO
            short_time_UT=change_direction_self_ind(find(time_after_CO_from_self_CD<=UT_time_from_CO));
            %2. find distnace of the bat from CO point
            Distance_UT_from_CO_self= pos_self_x(change_direction_self_ind)-pos_self_x(CO_point_before_UT(UT_i));
            close_UT=change_direction_self_ind(find(Distance_UT_from_CO_self<UT_distance_from_CO));
            
            if ismember(short_time_UT,close_UT)
                ii=ii+1;
                UT_point(ii)=short_time_UT(ismember(short_time_UT(1),close_UT));
                %                 UT_ind_in_CO=[UT_ind_in_CO CO_i];
                % find the start:
                long_dist_ind=find(abs(distnace_other_from_self)>=CO_window(1));
                before_CO=long_dist_ind(find([long_dist_ind-CO_point_before_UT(UT_i)]<0));
                start_ind=before_CO(end);%first ind before CO that the distance between the bats is lrager than ..
                if (CO_point_before_UT(UT_i)-start_ind)>max_wind_CO
                    start_ind=CO_point_before_UT(UT_i)-max_wind_CO;
                end
                UT_event_self{ii}=start_ind:UT_point(ii)+UT_window;
                UT_ind=[UT_ind UT_event_self{ii}];
            end
            
            %UT other:
            %1. find where in time there were UT
            time_after_CO_from_other_CD=change_direction_other_ind-CO_point_before_UT(UT_i);
            time_after_CO_from_other_CD(time_after_CO_from_other_CD<0)=NaN; %direction changed before CO
            short_time_UT=change_direction_other_ind(find(time_after_CO_from_other_CD<=UT_time_from_CO));
            %2. find distnace of the bat from CO point
            Distance_UT_from_CO_other= pos_other_x(change_direction_other_ind)-pos_other_x(CO_point_before_UT(UT_i));
            close_UT=change_direction_other_ind(find(Distance_UT_from_CO_other<UT_distance_from_CO));
            
            if ismember(short_time_UT,close_UT)
                ii=ii+1;
                UT_point(ii)=short_time_UT(ismember(short_time_UT(1),close_UT));
                %UT_ind_in_CO=[UT_ind_in_CO CO_i];
                % find the start:
                long_dist_ind=find(abs(distnace_other_from_self)>=CO_window(1));
                before_CO=long_dist_ind(find([long_dist_ind-CO_point_before_UT(UT_i)]<0));
                start_ind=before_CO(end);%first ind before CO that the distance between the bats is lrager than ..
                if (CO_point_before_UT(UT_i)-start_ind)>max_wind_CO
                    start_ind=CO_point_before_UT(UT_i)-max_wind_CO;
                end
                UT_event_other{ii}=start_ind:UT_point(ii)+UT_window;
                UT_ind=[UT_ind UT_event_other{ii}];
            end
        end
    end
end

% % Remove UT from CO:
% CO_include_UT=CO_ind;
% ind_to_remove=find(ismember(CO_ind,UT_ind));
CO_no_UT=CO_ind;
% CO_no_UT(ind_to_remove)=[];
%
%behavioral_modes.CO_include_UT=CO_include_UT;
behavioral_modes.CO_no_UT=CO_no_UT;
behavioral_modes.UT_ind=UT_ind;

behav=4;
behav_modes_plot(behav).name='CO';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(CO_no_UT,length(pos_self_x));
behav=5;
behav_modes_plot(behav).name='UT';
[~,behav_modes_plot(behav).start,behav_modes_plot(behav).end]=find_length_of_consecutive_ind(UT_ind,length(pos_self_x));

% %% Find obstacle ind:
% %1. find where the bat is closest to the obstacle in some sort distance
% close_ind=find(abs(pos_self_x-obstacle_pos)< close_to_obstacle);
% %2. find only the relevant crosses (session 2):
% relevant_time_ind=find((ts/1e6)>=obstacle_time);
% close_ind(find(close_ind<relevant_time_ind(1)))=[];
% 
% %3. find the consecutive ind:
% [length_obs,start_obs,end_obs]=find_length_of_consecutive_ind(close_ind,length(pos_self_x));
% obstacle_point=[];
% 
% %4. find the closest point to the obstacle:
% for obs_i=1:length(length_obs)
%     relevant_ind=start_obs(obs_i):end_obs(obs_i);
%     [~, closest_ind]=min([abs(pos_self_x(relevant_ind)-obstacle_pos)]);
%     obstacle_point=[obstacle_point relevant_ind(closest_ind)];
% end
% behavioral_modes.obstacle_point=obstacle_point;
% 
% % for control find obstacle point in session 1:
% %----------------------------------------------------------
% %1. find where the bat is closest to the obstacle in some sort distance
% close_ind=find(abs(pos_self_x-obstacle_pos)< close_to_obstacle);
% %2. find only the relevant crosses (session 1):
% relevant_time_ind=find((ts/1e6)<=obstacle_time);
% close_ind(find(close_ind<relevant_time_ind(1)))=[];
% 
% %3. find the consecutive ind:
% [length_obs,start_obs,end_obs]=find_length_of_consecutive_ind(close_ind,length(pos_self_x));
% obstacle_point=[];
% 
% %4. find the closest point to the obstacle:
% for obs_i=1:length(length_obs)
%     relevant_ind=start_obs(obs_i):end_obs(obs_i);
%     [~, closest_ind]=min([abs(pos_self_x(relevant_ind)-obstacle_pos)]);
%     obstacle_point=[obstacle_point relevant_ind(closest_ind)];
% end
% behavioral_modes.obstacle_point_session1=obstacle_point;
% 
% %compute velocity in 5 last flight of session 1
% mat= repmat([-300:300]',1,length(5));
% obs_point_idx = repmat(obstacle_point(end-5:end-1),[601,1]);
% obs_trigger_idx = obs_point_idx+mat;
% if max(obs_trigger_idx(:))>length(velocity_self_xy)
%     obs_trigger_idx(:,end)=[];
% end
% vel_trig_obs=nanmean(abs(velocity_self_xy(obs_trigger_idx)),2);
% behavioral_modes.vel_trig_no_obs_5_last=vel_trig_obs;
%% PLOT
figure('units','normalized','outerposition',[0 0 1 1])

% if tag_i==2
color_self=[0 1 0];
color_other=[1 0 0];

% title:
axes('position', [0.4, y_position(1), 0.1,  0.1]);
title(sprintf('Bat %d day %d - Total recording time %.1f min, Total flight time %.1f min',bat,day,(total_time/us_factor)/60,(length(FE_ind)/frame_per_second)/60), 'FontSize', 11, 'FontWeight', 'Bold');
box off;
axis off;
% else
%     color_self=[1 0 0];
%     color_other=[0 1 0];
% end
% All behavior
%--------------------------------------------------
behav_color=[0 0 1;0 1 0;1 0 0;0.5 0.5 0.5;1 1 0];

for time_i=1:length(time_vec_for_plot)-1
    ax(time_i)=axes('position',[x_position y_position(time_i) panel_size]);
    relevant_behav_ind=find(ts>=time_vec_for_plot(time_i) & ts<time_vec_for_plot(time_i+1));
    by_pass=behavioral_modes.bypass_ind(find(ts(behavioral_modes.bypass_ind)>=time_vec_for_plot(time_i) & ts(behavioral_modes.bypass_ind)<time_vec_for_plot(time_i+1)));
    plot(ts(relevant_behav_ind),pos_other_x(relevant_behav_ind),'.','color',color_other)
    hold on;
    plot(ts(relevant_behav_ind),pos_self_x(relevant_behav_ind),'.','color',color_self);
    hold on;
    if sum(isnan(pos_self_x(relevant_behav_ind)))==length(pos_self_x(relevant_behav_ind))
       set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))]) 
    else
    set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))],'ylim' ,[min(pos_self_x(relevant_behav_ind)) max(pos_self_x(relevant_behav_ind))])
    end
    set(gca,'XTick',[min(ts) max(ts)],'XTickLabel',round(([min(ts) max(ts)])/1e3/60))
    if ~isempty(by_pass)
        plot(ts(by_pass),pos_self_x(by_pass),'*k');
    end
    
    for behav_mod_i=1:length(behav_modes_plot)
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
        end
    end
end

%% manual corrections
if correct_manually
    % 1. Add CO
    q1=str2num(cell2mat(inputdlg('Do you want to add CO? (1/0)')));
    if q1
        how_many_CO=str2num(cell2mat(inputdlg('How many CO do you want to add?')));
        [co_time_to_add,~]=ginput(how_many_CO);
        [~,ind]=min(abs(ts-co_time_to_add'));
        [val,indx]=min(abs(distance_change_sign-ind));
        co_ind_to_add=distance_change_sign(indx(val<100)); %add CO only if it found closest point of changed sign
        behavioral_modes.CO_point=sort([CO_point co_ind_to_add']);
        for co_i=1:length(co_ind_to_add)
            co=co_ind_to_add(co_i);
            %find ax:
            poss_ax=find([time_vec_for_plot]<ts(co));
            axs=poss_ax(end);
            relevant_ax=ax(axs);
            axes(relevant_ax)
            plot(ts(co),pos_self_x(co),'ko')
            
        end
        manually_added.co_ind_to_add=co_ind_to_add;
    else
        manually_added.co_ind_to_add=[];
    end
    
    %2. Add solo
    q2=str2num(cell2mat(inputdlg('Do you want to correct solo? (1/0)')));
    if q2
        how_many_solo=str2num(cell2mat(inputdlg('How many Solo do you want to correct?')));
        manually_added.solo_ind_to_add=[];
        for solo_i=1:how_many_solo
            [solo_time_to_add,~]=ginput(2);
            [~,ind]=min(abs(ts-solo_time_to_add'));
            relevant_ind_to_add=ind(1):ind(2);
            relevant_ind_to_add=setdiff(relevant_ind_to_add,behavioral_modes.solo_ind);
            behavioral_modes.solo_ind=sort([behavioral_modes.solo_ind relevant_ind_to_add]);
            manually_added.solo_ind_to_add=[manually_added.solo_ind_to_add,relevant_ind_to_add];
            %find ax:
            poss_ax=find([time_vec_for_plot]<relevant_ind_to_add(1));
            axs=poss_ax(end);
            relevant_ax=ax(axs);
            axes(relevant_ax)
            plot([ts(relevant_ind_to_add(1)) ts(relevant_ind_to_add(end))],[ pos_self_x(relevant_ind_to_add(1)) pos_self_x(relevant_ind_to_add(end))],'k+')
        end
    else
        manually_added.solo_ind_to_add=[];
    end
    %2. Add tracking
    q3=str2num(cell2mat(inputdlg('Do you want to correct tracking? (1/0)')));
    if q3
        how_many_tr=str2num(cell2mat(inputdlg('How many tracking do you want to correct?')));
        manually_added.tracking_ind_to_add=[];
        
        for solo_i=1:how_many_tr
            [tr_time_to_add,~]=ginput(2);
            [~,ind]=min(abs(ts-tr_time_to_add'));
            relevant_ind_to_add=ind(1):ind(2);
            relevant_ind_to_add=setdiff(relevant_ind_to_add,behavioral_modes.tracking_ind);
            manually_added.tracking_ind_to_add=[manually_added.tracking_ind_to_add,relevant_ind_to_add];
            
            %find ax:
            poss_ax=find([time_vec_for_plot]<ts(relevant_ind_to_add(1)));
            axs=poss_ax(end);
            relevant_ax=ax(axs);
            axes(relevant_ax)
            plot([ts(relevant_ind_to_add(1)) ts(relevant_ind_to_add(end))],[ pos_self_x(relevant_ind_to_add(1)) pos_self_x(relevant_ind_to_add(end))],'k+')
            behavioral_modes.tracking_ind=sort([behavioral_modes.tracking_ind relevant_ind_to_add]);
            
        end
    else
        manually_added.tracking_ind_to_add=[];
    end
    file_name=fullfile(behave_struct_folder,['manually_added_behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
    save(file_name,'manually_added')
else
    file_name=fullfile(behave_struct_folder,['manually_added_behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
    load(file_name)
    %add:
    relev_fields=fields(manually_added);
    if strcmp(relev_fields,'tracking_ind_to_add')
        behavioral_modes.tracking_ind=sort([behavioral_modes.tracking_ind manually_added.tracking_ind_to_add]);
    end
    if strcmp(relev_fields,'solo_ind_to_add')
        behavioral_modes.solo_ind=sort([behavioral_modes.tracking_ind manually_added.solo_ind_to_add]);
    end
    if strcmp(relev_fields,'co_ind_to_add')
        behavioral_modes.CO_point=sort([behavioral_modes.CO_point manually_added.co_ind_to_add]);
    end
end
%% create time scale:
axes('position',[x_position 0.03 panel_size(1) 0.02]);
time_bin_in_min=(time_bin/us_factor)/60;
plot([0 time_bin/time_bin_in_min], [0 0],'k','LineWidth',5)
xlim([0 time_bin])
text(time_bin/(2*time_bin_in_min),-1.5, '1 min','HorizontalAlignment','center')
box off
axis off
%create legend
axes('position',[0.7 0.95 0.2 0.05])
x_pos=0:1/length(behav_modes_plot):1;
rec_size=0.8*(1/length(behav_modes_plot));
for behav_mod_i=1:length(behav_modes_plot)
    rec_x=[x_pos(behav_mod_i) x_pos(behav_mod_i)+rec_size x_pos(behav_mod_i)+rec_size x_pos(behav_mod_i)];
    rec_y=[0  0 1 1];
    p=patch(rec_x,rec_y,behav_color(behav_mod_i,:),'EdgeColor','none');
    set(p,'FaceAlpha',0.3)
    if behav_mod_i==3
        text(rec_x(1)+0.01,0.7,behav_modes_plot(behav_mod_i).name(1:5))
        text(rec_x(1)+0.01,0.3,behav_modes_plot(behav_mod_i).name(7:13))
    else
        text(rec_x(1)+0.01,0.5,behav_modes_plot(behav_mod_i).name)
    end
    box off
    axis off
end

% save
fig_name=fullfile(behave_analysis_fig_dir_out,['behavioral_mode_bat_',num2str(bat),'_day_',num2str(day),'.tif']);
saveas(gcf,fig_name)
fig_name=fullfile(behave_analysis_fig_dir_out,['behavioral_mode_bat_',num2str(bat),'_day_',num2str(day),'.fig']);
saveas(gcf,fig_name)
clf

%% Plot 2: behavioral measurments:
% figure param:
panel_size=[0.1 0.085];

%x pos
horizontal_dist=0.15;
x_position(1)=0.05;
x_position(2)=x_position(1)+horizontal_dist;
x_position(3)=x_position(2)+horizontal_dist;
x_position(4)=x_position(3)+horizontal_dist;
x_position(5)=x_position(4)+horizontal_dist;

%y pos
vertical_dist=0.15;
y_position(1)=0.84;
y_position(2)=y_position(1)-vertical_dist;
y_position(3)=y_position(2)-vertical_dist;
y_position(4)=y_position(3)-vertical_dist;
y_position(5)=y_position(4)-vertical_dist;
y_position(6)=y_position(5)-vertical_dist;
y_position(7)=y_position(6)-vertical_dist;
y_position(8)=y_position(7)-vertical_dist;
y_position(9)=y_position(8)-vertical_dist;

% space bins:
X_min=3;
X_max=120;
X_bin_size=5;
X_bins_vector=X_min:X_bin_size:X_max;
X_bins_vector_of_centers=X_bins_vector(1:end-1)+X_bin_size/2;
% space bins CO:
X_bin_size_CO=10;
X_bins_vector_CO=X_min:X_bin_size_CO:X_max;
X_bins_vector_of_centers_CO=X_bins_vector_CO(1:end-1)+X_bin_size_CO/2;

%Time bins:
T_min=ts(FE_ind(1));
T_max=ts(FE_ind(end));
T_bin_size=total_time/6; %
T_bins_vector=T_min:T_bin_size:T_max;
T_bins_vector_of_centers=T_bins_vector(1:end-1)+T_bin_size/2;

%% plot
% title:
axes('position', [0.4, y_position(1)+0.01, 0.1,  0.1]);
title(sprintf('Bat %d day %d - Total recording time %.1f min, Total flight time %.1f min',bat,day,(total_time/us_factor)/60,(length(FE_ind)/frame_per_second)/60), 'FontSize', 11, 'FontWeight', 'Bold');
box off;
axis off;


%1. Bar graph of time spent in different behaviors
%---------------------------------------------------------------------

% Precentage of flight time
ax1=axes('position',[x_position(1) y_position(1) panel_size]);
behavior_ind=[length(solo_ind),length(tracking_ind),length(CO_no_UT)];
h=100*(behavior_ind./length(FE_ind));
h(end+1)=100-sum(h);
bar(h)
set(gca,'XTick',1:length(behavior_ind)+1,'XTickLabel',{'Solo','Tracking','CO','None'})
ylabel('% of flight time')
title('Time spent in different behaviors')
xtickangle(ax1,45)
ylim([0 max(h)*1.1])
behavioral_modes.time_spent_in_different_behaviors=h;
% time histograms:
axes('position',[x_position(1) y_position(2) panel_size]);
h=hist(ts(solo_ind),T_bins_vector_of_centers);
sharp_ind_solo_hist_time=max(h)./mean(h);
%h=h/frame_per_second;
bar(T_bins_vector_of_centers-ts(FE_ind(1)),h)
set(gca,'XTick',[T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]);
set(gca,'XTickLabel',round([T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]/us_factor/60));
title(sprintf('Solo timings SI=%.2f',sharp_ind_solo_hist_time))
xlabel('Time (minute)')
ylabel('Count')
behavioral_modes.solo_hist_time=h;

axes('position',[x_position(1) y_position(3) panel_size]);
h=hist(ts(tracking_ind),T_bins_vector_of_centers);
sharp_ind_tracking_hist_time=max(h)./mean(h);
bar(T_bins_vector_of_centers-ts(FE_ind(1)),h)
set(gca,'XTick',[T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]);
set(gca,'XTickLabel',round([T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]/us_factor/60));
title(sprintf('tracking timings SI=%.2f',sharp_ind_tracking_hist_time))
xlabel('Time (minute)')
ylabel('Count')
behavioral_modes.tracking_hist_time=h;

axes('position',[x_position(1) y_position(5) panel_size]);
h=hist(ts(CO_point),T_bins_vector_of_centers);
sharp_ind_CO_hist_time=max(h)./mean(h);
%h=h/frame_per_second;
bar(T_bins_vector_of_centers-ts(FE_ind(1)),h)
set(gca,'XTick',[T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]);
set(gca,'XTickLabel',round([T_bins_vector_of_centers(1):T_bin_size:T_bins_vector_of_centers(end)]/us_factor/60));
title(sprintf('CO timings SI=%.2f',sharp_ind_CO_hist_time))
xlabel('Time (minute)')
ylabel('Count')
behavioral_modes.CO_hist_time=h;

%%
%2. histogrmas of time spent in different behaviors along the tunnle:
%---------------------------------------------------------------------

%a. solo
axes('position',[x_position(2) y_position(2) panel_size]);
h=hist(pos_self_x(solo_ind),X_bins_vector_of_centers);
h=h/frame_per_second;
sharp_ind_solo_hist_space=max(h)./mean(h);
bar(X_bins_vector_of_centers,h)
set(gca,'Xlim',[X_min X_max]);
title(sprintf('Solo position SI=%.2f',sharp_ind_solo_hist_space))
xlabel('X position (m)')
ylabel('Time (s)')
behavioral_modes.solo_hist_space=h;

%b. tracking
axes('position',[x_position(2) y_position(3) panel_size]);
h=hist(pos_self_x(tracking_ind),X_bins_vector_of_centers);
h=h/frame_per_second;
sharp_ind_tracking_hist_space=max(h)./mean(h);
bar(X_bins_vector_of_centers,h)
set(gca,'Xlim',[X_min X_max]);
title(sprintf('Tracking position SI=%.2f',sharp_ind_tracking_hist_space))
xlabel('X position (m)')
ylabel('Time (s)')
behavioral_modes.tracking_hist_space=h;

%c. tracked
% axes('position',[x_position(2) y_position(4) panel_size]);
% h=hist(pos_self_x(being_tracked_ind),X_bins_vector_of_centers);
% h=h/frame_per_second;
% bar(X_bins_vector_of_centers,h)
% set(gca,'Xlim',[X_min X_max]);
% title('Being tracked position')
% xlabel('X position (m)')
% ylabel('Time (s)')

%d. CO
axes('position',[x_position(2) y_position(5) panel_size]);
h=hist(pos_self_x(CO_point),X_bins_vector_of_centers_CO);
h_no_edges=h(1+bins_to_remove_from_edge_CO_hist:end-bins_to_remove_from_edge_CO_hist);
sharp_ind_CO_hist_space=max(h_no_edges)./mean(h_no_edges);
bar(X_bins_vector_of_centers_CO,h)
set(gca,'Xlim',[X_min X_max]);
title(sprintf('CO position SI=%.2f',sharp_ind_CO_hist_space))
xlabel('X position (m)')
ylabel('Count')
behavioral_modes.CO_hist_space=h;
%e. UT
% axes('position',[x_position(2) y_position(6) panel_size]);
% h=hist(pos_self_x(UT_point),X_bins_vector_of_centers);
% bar(X_bins_vector_of_centers,h)
% set(gca,'Xlim',[X_min X_max]);
% title('UT position')
% xlabel('X position (m)')
% ylabel('Count')

% %3. velocity X:
% %---------------------------------------------------------------------
% %a. velocity profile along the tunnle:
% axes('position',[x_position(2) y_position(1) panel_size]);
% plot(pos_self_x, abs(velocity_self),'.')
% set(gca,'Xlim',[X_min X_max]);
%
% %b. velocity triggered by obstacle
%
% %c. bar of mean velocity during: ALL solo tracked tracking
% ax2=axes('position',[x_position(2) y_position(2) panel_size]);
% mean_vel_behavior=[nanmean(abs(velocity_self(FE_ind))),nanmean(abs(velocity_self(solo_ind))),nanmean(abs(velocity_self(tracking_other_ind))),nanmean(abs(velocity_self(being_tracked_ind)))];
% [sem(1)]=fn_compute_sem(abs(velocity_self(FE_ind)));
% [sem(2)]=fn_compute_sem(abs(velocity_self(solo_ind)));
% [sem(3)]=fn_compute_sem(abs(velocity_self(tracking_other_ind)));
% [sem(4)]=fn_compute_sem(abs(velocity_self(being_tracked_ind)));
% bar(1:length(mean_vel_behavior),mean_vel_behavior); hold on;
% errorbar(1:length(mean_vel_behavior),mean_vel_behavior,sem,'.')
% set(gca,'XTick',1:4,'XTickLabel',{'All flight','Solo','Tracking','Tracked'})
% xtickangle(ax2,45)
% title('X Speed (m/s) in differnt behaviors')
%
% %d. velocity triggered by CO
% axes('position',[x_position(2) y_position(5) panel_size]);
% mat= repmat([-300:300]',1,length(CO_point));
% CO_point_idx = repmat(CO_point,[601,1]);
% CO_trigger_idx = CO_point_idx+mat;
% vel_trig_CO=nanmean(abs(velocity_self(CO_trigger_idx)),2);
% plot(-300:300,vel_trig_CO,'k')
% hold on;
% plot([0 0],[min(vel_trig_CO) max(vel_trig_CO)],'r')
% ylim([min(vel_trig_CO) max(vel_trig_CO)])
% xlim([-300 300])
% labels=get(gca,'Xtick');
% set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
% xlabel('Time (s)')
% ylabel('X Speed (m/s)')
% title('mean velocity triggered by CO')
% %d. velocity triggered by UT
% axes('position',[x_position(2) y_position(6) panel_size]);
% mat= repmat([-300:300]',1,length(UT_point));
% UT_point_idx = repmat(UT_point,[601,1]);
% UT_trigger_idx = UT_point_idx+mat;
% vel_trig_UT=nanmean(abs(velocity_self(UT_trigger_idx)),2);
% plot(-300:300,vel_trig_UT,'k')
% hold on;
% plot([0 0],[min(vel_trig_UT) max(vel_trig_UT)],'r')
% ylim([min(vel_trig_UT) max(vel_trig_UT)])
% xlim([-300 300])
% labels=get(gca,'Xtick');
% set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
% xlabel('Time (s)')
% ylabel('X Speed (m/s)')
% title('mean velocity triggered by UT')

%3. velocity XY:
%---------------------------------------------------------------------
%a. velocity profile along the tunnle:
axes('position',[x_position(3) y_position(1) panel_size]);
plot(pos_self_x, velocity_self_xy,'.')
set(gca,'Xlim',[X_min X_max]);
set(gca,'Ylim',[0 15]);
title('Velocity along the tunnle')
xlabel('X position (m)')
ylabel('Velocity (m/s)')
%b. velocity triggered by obstacle

%c. bar of mean velocity during: ALL solo tracked tracking
ax2=axes('position',[x_position(3) y_position(2) panel_size]);
%mean_vel_behavior=[nanmean(abs(velocity_self_xy(FE_ind))),nanmean(abs(velocity_self_xy(solo_ind))),nanmean(abs(velocity_self_xy(tracking_other_ind))),nanmean(abs(velocity_self_xy(being_tracked_ind)))];
mean_vel_behavior=[nanmean(abs(velocity_self_xy(FE_ind))),nanmean(abs(velocity_self_xy(solo_ind))),nanmean(abs(velocity_self_xy(tracking_ind)))];
[sem_vel(1)]=fn_compute_sem(abs(velocity_self_xy(FE_ind)));
[sem_vel(2)]=fn_compute_sem(abs(velocity_self_xy(solo_ind)));
[sem_vel(3)]=fn_compute_sem(abs(velocity_self_xy(tracking_ind)));
%[sem(4)]=fn_compute_sem(abs(velocity_self_xy(being_tracked_ind)));
bar(1:length(mean_vel_behavior),mean_vel_behavior); hold on;
errorbar(1:length(mean_vel_behavior),mean_vel_behavior,sem_vel,'.')
set(gca,'XTick',1:3,'XTickLabel',{'All flight','Solo','Tracking'})
xtickangle(ax2,45)
title('velocity (m/s) in differnt behaviors')
behavioral_modes.sem_vel=sem_vel;
behavioral_modes.mean_vel_behavior=mean_vel_behavior;

%d. velocity triggered by CO
axes('position',[x_position(3) y_position(5) panel_size]);
mat= repmat([-300:300]',1,length(CO_point));
CO_point_idx = repmat(CO_point,[601,1]);
CO_trigger_idx = CO_point_idx+mat;
if max(max(CO_trigger_idx))>length(velocity_self_xy)
  CO_trigger_idx(:,end)=[];  
end
vel_trig_CO=nanmean(abs(velocity_self_xy(CO_trigger_idx)),2);
plot(-300:300,vel_trig_CO,'k')
hold on;
plot([0 0],[min(vel_trig_CO) max(vel_trig_CO)],'r')
ylim([min(vel_trig_CO) max(vel_trig_CO)])
xlim([-300 300])
%ylim([6.5 8.5])
labels=get(gca,'Xtick');
set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
xlabel('Time (s)')
ylabel('X Speed (m/s)')
title('mean velocity triggered by CO')
behavioral_modes.vel_trig_CO=vel_trig_CO;
% 
% %d. velocity triggered by obstacle
% obs_point=behavioral_modes.obstacle_point;
% axes('position',[x_position(3) y_position(6) panel_size]);
% mat= repmat([-300:300]',1,length(obs_point));
% obs_point_idx = repmat(obs_point,[601,1]);
% obs_trigger_idx = obs_point_idx+mat;
% if max(obs_trigger_idx(:))>length(velocity_self_xy)
%     obs_trigger_idx(:,end)=[];
% end
% vel_trig_obs=nanmean(abs(velocity_self_xy(obs_trigger_idx)),2);
% plot(-300:300,vel_trig_obs,'k')
% hold on;
% plot([0 0],[min(vel_trig_obs) max(vel_trig_obs)],'r')
% ylim([min(vel_trig_obs) max(vel_trig_obs)])
% xlim([-300 300])
% labels=get(gca,'Xtick');
% set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
% xlabel('Time (s)')
% ylabel('X Speed (m/s)')
% title('mean velocity triggered by obstacle')
% behavioral_modes.vel_trig_obs=vel_trig_obs;
% behavioral_modes.vel_trig_obs_5_first=nanmean(abs(velocity_self_xy(obs_trigger_idx(:,1:5))),2);


% %d. velocity triggered by UT
% axes('position',[x_position(2) y_position(6) panel_size]);
% mat= repmat([-300:300]',1,length(UT_point));
% UT_point_idx = repmat(UT_point,[601,1]);
% UT_trigger_idx = UT_point_idx+mat;
% vel_trig_UT=nanmean(abs(velocity_self_xy(UT_trigger_idx)),2);
% plot(-300:300,vel_trig_UT,'k')
% hold on;
% plot([0 0],[min(vel_trig_UT) max(vel_trig_UT)],'r')
% ylim([min(vel_trig_UT) max(vel_trig_UT)])
% xlim([-300 300])
% labels=get(gca,'Xtick');
% set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
% xlabel('Time (s)')
% ylabel('X Speed (m/s)')
% title('mean velocity triggered by UT')



%4. Y position
%---------------------------------------------------------------------
%a. XY time spent
axes('position',[x_position(4) y_position(1) panel_size]);
y_bin=-1.5:0.25:1.5;
h=hist(pos_self_y,y_bin);
h=100*(h./sum(h));
bar(y_bin,h)
xlim([min(y_bin) max(y_bin)])
ylim([0 max(h)])
xlabel('Y position (m)')
ylabel('% of time')
title('Y time spent')

% axes('position',[x_position(3) y_position(1) panel_size]);
% plot(pos_self_x,pos_self_y,'.')
% set(gca,'Xlim',[X_min X_max]);
% set(gca,'ylim',[min(pos_self_y) max(pos_self_y)]);
% title('XY position')
% xlabel('X position (m)')
% ylabel('Y position (m)')

%b. distance in y during tracking
axes('position',[x_position(4) y_position(3) panel_size]);
y_dist_bin=-2:0.25:2;
h=hist(distnace_other_from_self_y(tracking_other_ind),y_dist_bin);
h=100*(h./sum(h));
bar(y_dist_bin,h)
xlim([min(y_dist_bin) max(y_dist_bin)])
xlabel('Y distance (m)')
ylabel('% of time')
title('Distance between the bats during tracking')

%c. distance in y during tracking
axes('position',[x_position(4) y_position(4) panel_size]);
y_dist_bin=-2:0.25:2;
h=hist(distnace_other_from_self_y(being_tracked_ind),y_dist_bin);
h=100*(h./sum(h));
bar(y_dist_bin,h)
xlim([min(y_dist_bin) max(y_dist_bin)])
xlabel('Y distance (m)')
ylabel('% of time')
title('Distance between the bats during being tracked')

%e. plot y pos triggered by CO
axes('position',[x_position(4) y_position(5) panel_size]);
for dir_i=1:2
    co_point_dir=intersect(CO_point,behavioral_modes.directional_ind{dir_i});
    mat= repmat([-300:300]',1,length(co_point_dir));
    CO_point_idx = repmat(co_point_dir',[601,1]);
    CO_trigger_idx = CO_point_idx+mat;
    if max(CO_trigger_idx(:))>length(velocity_self_xy)
        CO_trigger_idx(:,end)=[];
    end
    Y_trig_CO=pos_self_y(CO_trigger_idx);
    plot(-300:300,Y_trig_CO,'color',direc_color(dir_i,:))
    hold on;
end
plot([0 0],[min(pos_self_y) max(pos_self_y)],'r')
ylim([-1.7 1.7])
xlim([-300 300])
[x,y]=ds2nfu([150 300],[-1.3 -1.3]);
annotation('arrow',x,y,'Color',direc_color(1,:))
[x,y]=ds2nfu([300 150],[-1.5 -1.5]);
annotation('arrow',x,y,'Color',direc_color(2,:))
labels=get(gca,'Xtick');
set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
xlabel('Time (s)')
ylabel('y position (m)')
title('Y position triggered by CO')
% 
% %e. plot y pos triggered by obstacle
% axes('position',[x_position(4) y_position(6) panel_size]);
% for dir_i=1:2
%     obs_point_dir=intersect(obs_point,behavioral_modes.directional_ind{dir_i});
%     mat= repmat([-300:300]',1,length(obs_point_dir));
%     obs_point_idx = repmat(obs_point_dir',[601,1]);
%     obs_trigger_idx = obs_point_idx+mat;
%     if max(obs_trigger_idx(:))>length(velocity_self_xy)
%         obs_trigger_idx(:,end)=[];
%     end
%     
%     
%     Y_trig_obs=pos_self_y(obs_trigger_idx);
%     plot(-300:300,Y_trig_obs,'color',direc_color(dir_i,:))
%     hold on;
% end
% plot([0 0],[min(pos_self_y) max(pos_self_y)],'r')
% ylim([-1.5 1.5])
% xlim([-300 300])
% labels=get(gca,'Xtick');
% set(gca,'Xtick',labels,'XtickLabel',labels/frame_per_second)
% xlabel('Time (s)')
% ylabel('y position (m)')
% title('Y position triggered by obstacle')

% save
fig_name=fullfile(behave_analysis_fig_dir_out,['behavioral_analysis_bat_',num2str(bat),'_day_',num2str(day),'.tif']);
saveas(gcf,fig_name)

clf

%% save struct:

file_name=fullfile(behave_struct_folder,['behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
save(file_name,'behavioral_modes')

