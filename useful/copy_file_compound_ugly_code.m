%copy compund cells:
clear
out_dir_name='D:\Ayelet\2bat_proj\Analysis\new_code\figures\per_field_time_analysis\place_cells\compound';
origin_dir_name='D:\Ayelet\2bat_proj\Analysis\new_code\figures\per_field_time_analysis\place_cells\';
dir_info=dir(origin_dir_name);
dir_info_names={dir_info.name};
dir_info_names([dir_info.isdir])=[];
cell_ind=regexp(dir_info_names{1},'cell_');
cell_names=cellfun(@(a) a(cell_ind+5:cell_ind+8),dir_info_names,'UniformOutput',false);

cell_co_solo_initial_analysis_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis';
files = dir(cell_co_solo_initial_analysis_struct_folder);
cell_co_solo_initial_struct_names = {files.name};
cell_co_solo_initial_struct_names([files.isdir])=[];
%% for each cell
for ii_cell = 1:length(cell_co_solo_initial_struct_names)
    ii_cell
    load(fullfile(cell_co_solo_initial_analysis_struct_folder,cell_co_solo_initial_struct_names{ii_cell}));
    cell_num=cell_co_solo_initial_analysis.exp_data.cell_num ;
    compound=0;
    for ii_dir=1:2
        per_field=cell_co_solo_initial_analysis.co(ii_dir).per_field_href ;
        for ii_field=1:length(per_field)
            try
                if per_field(ii_field).time_signif_field.compound==1
                    compound=1;
                    Match=cellfun(@(a) find(contains(a,num2str(cell_num))) , cell_names, 'UniformOutput', 0);
                    r=find(cellfun(@(c) ~isempty(c),Match));
                    if ~isempty(r)
                    file_name=fullfile(origin_dir_name,dir_info_names{r});
                    copy_file_to_new_dir(file_name,out_dir_name)
                    end
                end
            catch
            end
            
        end
    end
    
end