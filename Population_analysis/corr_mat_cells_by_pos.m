function [mat_corr mean_corr]=corr_mat_cells_by_pos(mat_1,mat_2,dim)
if dim==2
    mat_1=mat_1';
    mat_2=mat_2';

end

corr_values=corr(mat_1,mat_2,'rows','pairwise');
mat_corr=diag(corr_values);
mean_corr=nanmean(mat_corr);