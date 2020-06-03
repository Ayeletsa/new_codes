function [cells_corr_between,cells_corr_within]=compute_corrs_within_and_between_cells(data)


cells_corr_within=[];
cells_corr_between=[];
for cell_i=1:length(data)
    
    mat_1=squeeze(data(cell_i,:,:));
    cell_corr_values=corr(mat_1(1,:)',mat_1(2,:)','rows','pairwise');
    cells_corr_within=[cells_corr_within,cell_corr_values];
    
    cell_vec=zeros(length(data),1);
    cell_vec(cell_i)=1;
    other_cells=data(~cell_vec,:,:);
    other_cells=[squeeze(other_cells(:,1,:));squeeze(other_cells(:,2,:))];
    other_cell_corr_values=corr(mat_1(1,:)',other_cells','rows','pairwise');
    cells_corr_between=[cells_corr_between;other_cell_corr_values(:)];
    other_cell_corr_values=corr(mat_1(2,:)',other_cells','rows','pairwise');
    cells_corr_between=[cells_corr_between;other_cell_corr_values(:)];
    
end
