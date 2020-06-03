function all_co_ind_to_add=mannually_correct_co(pos_self_x,ts,distance_change_sign,manual_min_dis_from_CO)


    all_co_ind_to_add=[];
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
            [co_time_to_add,~]=ginput(1);
            %find closeset point:
            [~,ind]=min(abs(ts-co_time_to_add'));
            [val,indx]=min(abs(distance_change_sign-ind));
            co_ind_to_add=distance_change_sign(indx(val<manual_min_dis_from_CO));
            %co_ind_to_add=ind;
            % plot:
            plot(ts(co_ind_to_add),pos_self_x(co_ind_to_add),'b*','markersize',20);

            %add to data:
             all_co_ind_to_add=[all_co_ind_to_add,co_ind_to_add];
        end
        % if 'd' --> delete last point and clear from figure:
        if strcmp(key_pressed,'d')
            all_co_ind_to_add(end)=[];
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