function bypass_ind=find_bypass_ind(distnace_other_from_self,tracking_ind)


% Defined by change from tracking to being tracked
%a. find where distance change signs:
distnace_other_from_self_shifted_plus=[0; distnace_other_from_self];
distnace_other_from_self_shifted_minus=[distnace_other_from_self; 0];
distance_change_sign=find(sign(distnace_other_from_self_shifted_minus.*distnace_other_from_self_shifted_plus)==-1)-1;
% find intersections with tracking:
bypass_ind=intersect(distance_change_sign,tracking_ind);