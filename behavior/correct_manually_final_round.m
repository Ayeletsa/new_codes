% code to correct co mannually final:
clear
behave_day_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\go_over_behavior_again\';
dir_info=dir(behave_day_struct_folder);
manual_min_dis_from_CO=100;
for day_i=18:length(dir_info)
    day_i
    %% load data:
    behav_file_name=fullfile(behave_day_struct_folder,dir_info(day_i).name);
    load(behav_file_name)
    %ugly!!
    bat=behav_file_name(129:132);
    day=behav_file_name(138:145);
    general_behavior_data_file_name=['D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\general_behave_analysis_bat',bat,'_day',day,'.mat\'];
    load(general_behavior_data_file_name)
    %% plot figure:
    behave_plot(behav_file_name,general_behavior_data_file_name)
    
    %% correct
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
            plot(ts(co_ind_to_add),pos_self_x(co_ind_to_add),'ko');
            plot(ts(co_ind_to_add),pos_self_x(co_ind_to_add),'*');

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
    
    
    %% add to data:
    behavioral_modes.CO_point=sort([behavioral_modes.CO_point all_co_ind_to_add]);    
    
    %% save data:
    %1. save all struct:
    save(behav_file_name,'behavioral_modes')
    %2. save just mannualy adda:
    manually_added_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\behavioral_modes\day_structs\manaully_added_final\';
    manually_added_file=[manually_added_folder,'manually_added_final_bat',bat,'_day',day,'.mat\'];
    save(manually_added_file,'all_co_ind_to_add')

    %3. save figure:
    manually_added_fig=[manually_added_folder,'manually_added_final_bat',bat,'_day',day,'.jpg'];
    saveas(gcf,manually_added_fig)

    % close figure
    close all
end

function key_switch(hObject, eventdata, handles)
KEY=eventdata.Key;
setappdata(hObject,'KEY',KEY);
end