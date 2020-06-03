function  rise_and_width_data=get_per_field_tuning_width_and_rise_time(per_field_time_signif_data,min_dis_pos_neg)

% get data of tuning width and rise time:
neg_width=per_field_time_signif_data.width_neg_interp;
neg_width(isnan(neg_width))=[];
pos_width=per_field_time_signif_data.width_pos_interp;
pos_width(isnan(pos_width))=[];
pos_rise_time=per_field_time_signif_data.pos_rise_time_interp;
pos_rise_time(isnan(pos_rise_time))=[];
neg_rise_time=per_field_time_signif_data.neg_rise_time_interp;
neg_rise_time(isnan(neg_rise_time))=[];

%define the tuning for coulring in plot if it is pos or neg:
pos_per_field_tuning=~isempty(pos_width);
neg_per_field_tuning=~isempty(neg_width);

%initialize for compound:
compound_width=[];
pos_compound=[];
neg_compound=[];

if ~isempty(pos_width) & ~isempty(neg_width) % if there are both positive and negative
    x_pos=per_field_time_signif_data.width_line_x_pos_interp;
    x_neg=per_field_time_signif_data.width_line_x_neg_interp;
    x_pos(isnan(x_pos))=[];x_neg(isnan(x_neg))=[];
    % check if there are clsoe to each
    % other:
    if abs(min(x_pos(:))-max(x_neg(:)))<=min_dis_pos_neg | abs(max(x_pos(:))-min(x_neg(:)))<=min_dis_pos_neg
        %TO DO: to correct it to work
        %properly on all type of tuning
        %curve if we will use it!
        compound_width=[max([x_pos(:);x_neg(:)])-min([x_pos(:);x_neg(:)])];
        pos_compound=sum(pos_width); %think how we want to compute this if there is more than one pos!
        neg_compound=sum(neg_width);
        pos_width=[];
        neg_width=[];
        neg_rise_time=[];
        pos_rise_time=[];
    end
end

%save to struct:
rise_and_width_data.pos_width=pos_width;
rise_and_width_data.neg_width=neg_width;
rise_and_width_data.compound_width=compound_width;
rise_and_width_data.pos_compound_width=pos_compound;
rise_and_width_data.neg_compound_width=neg_compound;
rise_and_width_data.neg_rise_time=neg_rise_time;
rise_and_width_data.pos_rise_time=pos_rise_time;
rise_and_width_data.pos_per_field_tuning=pos_per_field_tuning;
rise_and_width_data.neg_per_field_tuning=neg_per_field_tuning;

end


