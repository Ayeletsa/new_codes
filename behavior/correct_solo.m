function solo_ind_to_remove=correct_solo(solo_start,solo_end,manual_min_dis_from_CO,ts,pos_self_x)   
    solo_ind_to_remove=[];
    f=gcf;
    f.KeyPressFcn = @key_switch;
    
    W=1;
    while W 
        setappdata(f,'KEY','1');
        pause(0.001)
        key_pressed=getappdata(f,'KEY');
        % if space --> next figure
        if strcmp(key_pressed,'space')
            break
        end
        % if 'c' correct co:
        if strcmp(key_pressed,'c')
            [solo_ts_to_remove,~]=ginput(2);
            %find closeset point:
            [~,ind]=min(abs(ts-solo_ts_to_remove'));
            
            [val_start,indx_start]=min(abs(solo_start-ind(1)));
            [val_end,indx_end]=min(abs(solo_end-ind(2)));

            solo_start_to_remove=solo_start(indx_start(val_start<manual_min_dis_from_CO));
            solo_end_to_remove=solo_end(indx_end(val_end<manual_min_dis_from_CO));

            %co_ind_to_add=ind;
            % plot:
            plot(ts([solo_start_to_remove solo_end_to_remove]),pos_self_x([solo_start_to_remove solo_end_to_remove]),'b*','markersize',20);

            %add to data:
             solo_ind_to_remove=[solo_ind_to_remove,solo_start_to_remove:solo_end_to_remove];
        end
        % if 'd' --> delete last point and clear from figure:
        if strcmp(key_pressed,'d')
            solo_ind_to_remove(end)=[];
            % delete from figure:
            children = get(gca, 'children');
            delete(children(1));
        end
     end
    
   
function key_switch(hObject, eventdata, handles)
KEY=eventdata.Key;
setappdata(hObject,'KEY',KEY);
end 

end