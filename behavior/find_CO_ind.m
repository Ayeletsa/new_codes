function [CO_ind,CO_point,CO_event]=find_CO_ind(general_behavior_data_file_name,behav_params_file_name)
load(general_behavior_data_file_name)
load(behav_params_file_name)

%Cross over is defined when both bats are flying in differnt direction for t time before and t time after
%the distance change sign and te bats are in different directions:

%1) Both bats are flying:
both_bats_are_flying=intersect(FE_ind,bsp_proc_data(3-tag_i).flight_ind);
both_bats_are_flying_FE=find(ismember(FE_ind,both_bats_are_flying));
%2) bats with different flight directions:
different_flight_direc_ind_FE=find(sign(velocity_other_FE.*velocity_self_FE)==-1);
[different_direc_length,different_direc_start,different_direc_end]=find_length_of_consecutive_ind(FE_ind(different_flight_direc_ind_FE),length(pos_self_x));

%3) intercet prev conditions:
relevant_flights_FE=intersect(both_bats_are_flying_FE,different_flight_direc_ind_FE);
%relevant_flights_FE=intersect(relevant_flights_FE,short_distance_ind_FE);

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