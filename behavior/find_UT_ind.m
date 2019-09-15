function [CO_no_UT,UT_ind]=find_UT_ind(general_behavior_data_file_name,behav_params_file_name,CO_ind)


load(general_behavior_data_file_name)
load(behav_params_file_name)


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

%different directions:
different_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==-1);
[different_direc_length,different_direc_start,different_direc_end]=find_length_of_consecutive_ind(FE_ind(different_flight_direc_ind_FE),length(pos_self_x));


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


CO_no_UT=CO_ind;
