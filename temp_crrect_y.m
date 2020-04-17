% correct y in cell structs:
% open correct and save cell structs:
dir_cell_structs='D:\Ayelet\Data\Data_Nlg_Proc\yr_2018_bat_2389\cell_structs';
dir_info=dir(dir_cell_structs);
%run over all cells:
for cell_i=3:length(dir_info)
    %1. load cell struct:
    file_name=fullfile(dir_cell_structs,dir_info(cell_i).name);
    load(file_name)
    %2. remove y our of tunnel:
    outliers_pos_dist_from_midline=1.5;
    for tag_i=1:2
        cell_struct.bsp_data(tag_i).outliers_pos_dist_from_midline=outliers_pos_dist_from_midline;%update parameter
        %find ind to remove:
        ind_to_remove_from_bsp=(cell_struct.bsp_data(tag_i).pos_upsampled(2,:)<-outliers_pos_dist_from_midline | cell_struct.bsp_data(tag_i).pos_upsampled(2,:)>outliers_pos_dist_from_midline);
        %remove from pos and ts:
        cell_struct.bsp_data(tag_i).pos_upsampled=cell_struct.bsp_data(tag_i).pos_upsampled(:,~ind_to_remove_from_bsp);
        cell_struct.bsp_data(tag_i).ts_us_upsampled=cell_struct.bsp_data(tag_i).ts_us_upsampled(~ind_to_remove_from_bsp);
    end
    %3. save cell struct
    save(file_name,'cell_struct')
end