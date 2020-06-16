
function copy_file_to_new_dir(file_name,out_dir_name)
if isempty(dir(out_dir_name))
    mkdir(out_dir_name);
end
copyfile(file_name,out_dir_name);