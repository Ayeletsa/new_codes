function [field_center,field_size,field_height,field_edges,PSTH] = find_fields_PSTH_Basic_withShuffle(flights_struct,solo_param_file_name)
% flights struct include the following 1 dimensional fieldss: 
%                       pos             1xN
%                       ts_nlg_usec     1xN
%                       spike_pos       1xM
% bin_edges 1xL
% ref_height_for_width =    percentage from peak height (e.g 0.2)in which to
%                           compute the width of thw field
% ref_height_for_overlap =  percentage from peak height (e.g 0.5), in which
%                           to consider two overlapping fields.

load(solo_param_file_name)

field_center = [];
field_size = [];
field_height = [];
field_edges = [];
field_ind = [];
field_edges_ind = [];


bsp_pos = [flights_struct(:).pos]; 
ts = [flights_struct(:).ts_nlg_usec];
spike_pos = [flights_struct(:).spike_pos];
ts(isnan(ts)) = []; bsp_pos(isnan(bsp_pos)) = [];

[PSTH,~,~] = computePSTH(bsp_pos,ts,spike_pos,bin_edges);

bin_size = (bin_edges(end) - bin_edges(1))/(length(bin_edges)-1);
bin_centers = (bin_edges(1)+bin_size/2):bin_size:bin_edges(end);

PSTH_Nan2zero = PSTH; PSTH_Nan2zero(isnan(PSTH)) = 0; %changing nans to 
% zeros so that border peaks could be found. 
[pks,loc_ind] = findpeaks(PSTH_Nan2zero,'minPeakHeight',minPeakFR);
loc = bin_centers(loc_ind);

if isempty(pks)
    return
end

%% find field width and remove fields with not enough spikes
[field_size_bin, field_edges_bin] = find_field_width_edges(PSTH,pks,loc_ind,ref_height_for_width);

field_size = bin_size*field_size_bin; %units m
field_edges(1,:) = interp1(1:length(bin_centers),bin_centers,field_edges_bin(1,:));
field_edges(2,:) = interp1(1:length(bin_centers),bin_centers,field_edges_bin(2,:));

toRemove = [];
i = 1;
while i <= length(pks)
    
    num_spikes_inField = sum(and(spike_pos>field_edges(1,i),spike_pos<field_edges(2,i)));
    if num_spikes_inField <= min_spikes_perField
        toRemove = [toRemove i];
    end
    i = i+1;
end

field_edges(:,toRemove) = [];
field_size(toRemove) = [];
field_center = loc; field_center(toRemove) = [];
field_height = pks; field_height(toRemove) = [];
field_ind = loc_ind; field_ind(toRemove) = [];
field_edges_bin(:,toRemove) = [];
field_edges_ind(1,:) = floor(field_edges_bin(1,:));
field_edges_ind(2,:) = ceil(field_edges_bin(2,:));

if isempty(field_center)
    return
end

%% remove overlapping, lower fields:
if length(field_center)>1 %more than 1 field
    [field_height_sorted,ind_sorted] = sort( field_height,'descend');
    
    if ref_height_for_width == ref_height_for_overlap
        field_edges_overlap = field_edges;
    else
        [~, field_edges_overlap_bin] = find_field_width_edges(PSTH,field_height,field_ind,ref_height_for_overlap);
        field_edges_overlap(1,:) = interp1(1:length(bin_centers),bin_centers,field_edges_overlap_bin(1,:));
        field_edges_overlap(2,:) = interp1(1:length(bin_centers),bin_centers,field_edges_overlap_bin(2,:));
    end
    
    field_edges_sorted = field_edges_overlap(:,ind_sorted);
    field_center_sorted = field_center(ind_sorted);
    field_ind_sorted = field_ind(ind_sorted);
    toRemove = zeros([1,length(field_center)]);
    
    i = 2;
    while i <= length(field_height_sorted)
        
        edges = field_edges_sorted(:,i);
        %check if there are overlapping fields among the higher fields (1:i).
        overlapping_fields_ind1 = or(and(edges(1)>field_edges_sorted(1,1:i),edges(1)<field_edges_sorted(2,1:i)),...
            and(edges(2)>field_edges_sorted(1,1:i),edges(2)<field_edges_sorted(2,1:i)));
        overlapping_fields_ind2 = or(and(field_edges_sorted(1,1:i)>edges(1),field_edges_sorted(1,1:i)<edges(2)),...
            and(field_edges_sorted(2,1:i)>edges(1),field_edges_sorted(2,1:i)<edges(2)));
        overlapping_fields_ind = find(or(overlapping_fields_ind1,overlapping_fields_ind2));
        
        if isempty(overlapping_fields_ind) %if no overlapping fields, dont remove current field
            i = i+1;
            continue;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         limit_field_sizeR = false; limit_field_sizeL = false;
%         overlapping_fields_centers = field_center_sorted(overlapping_fields_ind);
%         temp = find(overlapping_fields_centers>field_center_sorted(i));
%         [~,ind_temp] = min(overlapping_fields_centers(temp));
%         overlapp_fieldR_ind = overlapping_fields_ind(temp(ind_temp));
%         height_overlap_field_R = field_height_sorted(overlapp_fieldR_ind);
%         
%         temp = find(overlapping_fields_centers<field_center_sorted(i));
%         [~,ind_temp] = max(overlapping_fields_centers(temp));
%         overlapp_fieldL_ind = overlapping_fields_ind(temp(ind_temp));
%         height_overlap_field_L = field_height_sorted(overlapp_fieldL_ind);
%         
%         if ~isempty(overlapp_fieldR_ind)
%             RelevantPSTH_R = (-1)*PSTH(field_ind_sorted(i):field_ind_sorted(overlapp_fieldR_ind)); %PSTH from first field peak to second field peak
%             [pks,loc_ind] = findpeaks(RelevantPSTH_R);
%             [valley_height_R,valley_ind_R] = min((-1)*pks);
%             valey_loc_R = bin_centers(field_ind_sorted(i) - 1 + loc_ind(valley_ind_R));
%             if valley_height_R<=(ref_height_for_width*height_overlap_field_R) %will go below 50% of the lower.....
%                 field_edges_sorted(2,i) = valey_loc_R;
%                 field_edges(2,ind_sorted(i)) = valey_loc_R;
%                 field_size(ind_sorted(i)) = diff(field_edges(:,ind_sorted(i)));
%                 limit_field_sizeR = true;
%             end
%         else
%             limit_field_sizeR = true;
%         end
%         
%         if ~isempty(overlapp_fieldL_ind)
%             RelevantPSTH_L = (-1)*PSTH(field_ind_sorted(overlapp_fieldL_ind):field_ind_sorted(i)); %PSTH from first field peak to second field peak
%             [pks,loc_ind] = findpeaks(RelevantPSTH_L);
%             [valley_height_L,valley_ind_L] = min((-1)*pks);
%             valey_loc_L = bin_centers(field_ind_sorted(overlapp_fieldL_ind) - 1 + loc_ind(valley_ind_L));
%             if valley_height_L<=(ref_height_for_width*height_overlap_field_L)
%                 field_edges_sorted(1,i) = valey_loc_L;
%                 field_edges(1,ind_sorted(i)) = valey_loc_L;
%                 field_size(ind_sorted(i)) = diff(field_edges(:,ind_sorted(i)));
%                 limit_field_sizeL = true;
%             end
%         else
%             limit_field_sizeL = true;
%         end
%   
%         if limit_field_sizeR && limit_field_sizeL
%             i = i+1;
%             continue;
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        toRemove(ind_sorted(i)) = 1;
        field_height_sorted(i) = [];
        field_edges_sorted(:,i) = [];
        field_center_sorted(i) = [];
        field_ind_sorted(i) = [];
        ind_sorted(i) = [];
    end
    
    field_height(toRemove>0) = [];
    field_size(toRemove>0) = [];
    field_center(toRemove>0) = [];
    field_edges(:,toRemove>0) = [];
    field_edges_ind(:,toRemove>0) = [];
    field_ind(toRemove>0) = [];
end

if isempty(field_center)
    return
end

%%
% % %% two peak firing rate threshold options:
% % % Option 1: the peak is above the 85% of firing rates of the PSTH - global 
% % % condition 
% % ind_cond1 = field_height >= TH_firingRate_percentiles;
% % % Option 2: the peak is above 1 Hz, and in a distance of 5 m from the 
% % % field's half height, limits the FR drops below 10% of the peak - local
% % % condition.
% % [~, field_edges_bin_05] = find_field_width_edges(PSTH,field_height,field_ind,0.5);
% % field_edges_bin_05(1,:) = floor(field_edges_bin_05(1,:));
% % field_edges_bin_05(2,:) = ceil(field_edges_bin_05(2,:));
% %     
% % dis = 5; %the condition is measured in a distance of 10 m from both field's edges.
% % dis_bin = dis/bin_size; %in bin units
% % Option2_TH = 0.1; %10% from field's peak
% % 
% % for i = 1:length(field_center)
% %     if field_edges_bin_05(1,i) > dis_bin
% %         PSTH_left = PSTH((field_edges_bin_05(1,i)-dis_bin):field_edges_bin_05(1,i));
% %     else
% %         PSTH_left = PSTH(1:field_edges_bin_05(1,i));
% %     end
% %     check_left = sum(PSTH_left<Option2_TH*field_height(i));
% %     
% %     if (field_edges_bin_05(2,i)+dis_bin)<=length(PSTH)
% %         PSTH_right = PSTH(field_edges_bin_05(2,i):(field_edges_bin_05(2,i)+dis_bin));
% %     else
% %         PSTH_right = PSTH(field_edges_bin_05(2,i):end);
% %     end
% %     check_right = sum(PSTH_right<Option2_TH*field_height(i));
% %     if (check_left>0)&&(check_right>0)
% %         ind_cond2(i) = true;
% %     else
% %         ind_cond2(i) = false;
% %     end
% % end
% % 
% % ind_cond_total = or(ind_cond1,ind_cond2);
% % toRemove = find(~ind_cond_total);
% % field_height(toRemove) = [];
% % field_size(toRemove) = [];
% % field_center(toRemove) = [];
% % field_edges(:,toRemove) = [];
% % field_edges_ind(:,toRemove) = [];
% % field_ind(toRemove) = [];
% % 
% % if isempty(field_center)
% %     return
% % end


%% limit field size if neighboring fields overlap (in field width and not 
%  in ref height for overlap)
[field_center,ind_sorted] = sort( field_center,'ascend');
field_height = field_height(ind_sorted);
field_size = field_size(ind_sorted);
field_ind = field_ind(ind_sorted);
field_edges = field_edges(:,ind_sorted);
field_edges_ind = field_edges_ind(:,ind_sorted);


if (ref_height_for_width < ref_height_for_overlap) &&(size(field_center,2)>1)
    checkEdgesOrder = (field_edges(1,2:end) - field_edges(2,1:end-1))<0;
    ind = find(checkEdgesOrder);
    %field i left border minus field i-1 right border. if it's  smalller
    %than 0 then, the field i overlap woth field i-1.
    if ~isempty(ind)
        for i=1:length(ind)
            fieldEdges1 = field_edges(:,ind(i));
            fieldEdges2 = field_edges(:,ind(i)+1);
            start_ind = field_ind(ind(i));
            end_ind = field_ind(ind(i)+1);
            RelevantPSTH = (-1)*PSTH(start_ind:end_ind);
            [pks,loc_ind] = findpeaks(RelevantPSTH);
            [~,max_pk_ind] = max(pks);
            valey_loc = bin_centers(start_ind - 1 + loc_ind(max_pk_ind));
            if valey_loc<fieldEdges1(2)
                fieldEdges1(2) = valey_loc;
            end
            if valey_loc>fieldEdges2(1)
                fieldEdges2(1) = valey_loc;
            end
            field_edges(:,ind(i)) = fieldEdges1;
            field_edges(:,ind(i) + 1) = fieldEdges2;
            field_size(ind(i)) = diff(fieldEdges1);
            field_size(ind(i) + 1) = diff(fieldEdges2);
        end
    end
end
        
        
        
%% remove unstable fields:

spikes_pos_per_flights = flights_struct;
num_of_flightsTH = max(num_of_flights_minimum,length(flights_struct)*stabilityTH);
toRemove = [];
i = 1;
while i<= length( field_center )
    count_flights = 0;
    edges = field_edges(:,i);
    
    for  j = 1:length(spikes_pos_per_flights)
        spike_pos_flight = spikes_pos_per_flights(j).spike_pos;
        isActiveInField = sum(and(spike_pos_flight >= edges(1), spike_pos_flight <= edges(2)))>0;
        count_flights = count_flights + isActiveInField;
        
        if count_flights >= num_of_flightsTH
            num_spikes_inField = sum(and(spike_pos>edges(1),spike_pos<edges(2)));
            if num_spikes_inField <= min_spikes_perField
                toRemove = [toRemove i];
            end
            i = i+1;
            break;
        end
    end
    if count_flights < num_of_flightsTH
        toRemove = [toRemove i];
        i = i+1;
    end
end

field_height(toRemove) = [];
field_size(toRemove) = [];
field_center(toRemove) = [];
field_edges(:,toRemove) = [];


%% Remove non significant fields based on shuffling analysis.
TH = 95;
signif_fields_ind = find_signifincant_fields(field_center, field_size, field_edges, flights_struct,TH);
field_center(~signif_fields_ind) = [];
field_size(~signif_fields_ind) = [];
field_height(~signif_fields_ind) = [];
field_edges(:,~signif_fields_ind) = [];
        
        

end