function [ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(vec_ind,length_long_vec)

%this function find the length of consecutive ind, returns:
%inputs: 
%vec_ind= indeces in the long vec of some event
%length_long_vec= length long vec of the entire time
%outpts:
%ind_length = vector with all consecutive ind length


    long_vec=zeros(1,length_long_vec);
    long_vec(vec_ind)=1;
    ind=diff(long_vec);
    start_ind=find(ind==1)+1;
    end_ind=find(ind==-1);
    % if tracking from first ind:
    if long_vec(1)==1
        start_ind=[1,start_ind];
    end
    
     if long_vec(end)==1
        end_ind=[end_ind length(long_vec)];
    end
    
    %find  length
    ind_length=end_ind-start_ind+1;
    