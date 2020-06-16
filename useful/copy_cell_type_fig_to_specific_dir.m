%% load dir info of figures:
origin_dir_name='D:\Ayelet\2bat_proj\Analysis\new_code\figures\basic_co_analysis\';
dir_info=dir(origin_dir_name);
dir_info_names={dir_info.name};
dir_info_names([dir_info.isdir])=[];
cell_ind=regexp(dir_info_names{1},'cell_');
cell_names=cellfun(@(a) a(cell_ind+5:cell_ind+8),dir_info_names,'UniformOutput',false);
%% load signif info:
inclusion_cells_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\inclusion_cells_struct';
filename = fullfile(inclusion_cells_folder_name,['all_cell_inclusion_list']);
load(filename)
%% transfer files:
%place cells files:
out_dir_name=[origin_dir_name,'\place_cells\'];

Match=cellfun(@(a) find(contains(a,strsplit(num2str(all_cell_inclusion_list.place_cells_signif)))) , cell_names, 'UniformOutput', 0);
r=find(cellfun(@(c) ~isempty(c),Match));
for ii_file=r
    file_name=fullfile(origin_dir_name,dir_info_names{ii_file});
    copy_file_to_new_dir(file_name,out_dir_name)
end

%ego cells files:
out_dir_name=[origin_dir_name,'\ego_cells\'];
Match=cellfun(@(a) find(contains(a,strsplit(num2str(all_cell_inclusion_list.ego_signif)))) , cell_names, 'UniformOutput', 0);
r=find(cellfun(@(c) ~isempty(c),Match));
for ii_file=r
    file_name=fullfile(origin_dir_name,dir_info_names{ii_file});
    copy_file_to_new_dir(file_name,out_dir_name)
end




