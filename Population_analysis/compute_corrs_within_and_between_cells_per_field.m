function [cells_corr_between,cells_corr_within,mean_corr_cell_value]=compute_corrs_within_and_between_cells_per_field(data)


cells_corr_within=[];
cells_corr_between=[];
for cell_i=1:length(data)
    cell_vec=zeros(length(data),1);
    cell_vec(cell_i)=1;
    mat_1=[data{cell_i,:}];
    cell_corr_values=corr(mat_1,mat_1,'rows','pairwise');
    cell_corr_values=triu(cell_corr_values,1);
    cell_corr_values=cell_corr_values(find(cell_corr_values~=0));
    cells_corr_within=[cells_corr_within; cell_corr_values];
    mean_corr_cell_value(cell_i)=nanmean(cell_corr_values);
    
    other_cells_mat=[data{~cell_vec,:}];
    
    other_cell_corr_values=corr(mat_1,other_cells_mat,'rows','pairwise');
    cells_corr_between=[cells_corr_between;other_cell_corr_values(:)];
end
