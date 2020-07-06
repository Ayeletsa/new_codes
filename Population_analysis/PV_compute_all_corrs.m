function corrs=PV_compute_all_corrs(ego_bin_start,PSTH_co_ego_bins,solo_PSTH,PSTH_solo_shuffle_ego_bins,solo_PSTH_other_dir,n_shuffles)

%1. average pos coding over cells:
% a. solo vs. co
for ego_bin_i=1:length(ego_bin_start)
    % a. solo vs. co
    co_allo_cell_mat=PSTH_co_ego_bins(:,:,ego_bin_i);
    solo_allo_cell_mat=solo_PSTH;
    [corrs.PV.co.mat(:,ego_bin_i), corrs.PV.co.tuning(ego_bin_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
    [corrs.cells.co.mat(:,ego_bin_i), corrs.cells.co.tuning(ego_bin_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,1);
    
    % b. solo vs. co norm max
    co_allo_cell_mat=co_allo_cell_mat./max(co_allo_cell_mat,[],2);
    solo_allo_cell_mat=solo_allo_cell_mat./max(solo_allo_cell_mat,2);
    [corrs.norm_PV.co.mat(:,ego_bin_i), corrs.norm_PV.co.tuning(ego_bin_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
    
    %c. other dir solo vs. co
    co_allo_cell_mat=PSTH_co_ego_bins(:,:,ego_bin_i);
    solo_allo_cell_mat=solo_PSTH_other_dir;
    [corrs.PV.other_dir.mat(:,ego_bin_i), corrs.PV.other_dir.tuning(ego_bin_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
    [corrs.cells.other_dir.mat(:,ego_bin_i), corrs.cells.other_dir.tuning(ego_bin_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,1);
    
    %d. other dir solo vs. co - norm to max
    co_allo_cell_mat=co_allo_cell_mat./max(co_allo_cell_mat,[],2);
    solo_allo_cell_mat=solo_allo_cell_mat./max(solo_allo_cell_mat,2);
    [corrs.norm_PV.other_dir.norm.mat(:,ego_bin_i), corrs.norm_PV.other_dir.tuning(ego_bin_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
    
    for shuffle_i=1:n_shuffles
        % Upper bound - solo (shuffle) vs. solo (all)
        %=====================================================
        % a. solo vs. co
        co_allo_cell_mat=PSTH_solo_shuffle_ego_bins(:,:,shuffle_i,ego_bin_i);
        solo_allo_cell_mat=solo_PSTH;
        [corrs.PV.co.shuffle.mat(:,ego_bin_i,shuffle_i), corrs.PV.co.shuffle.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
        [corrs.cells.co.shuffle.mat(:,ego_bin_i,shuffle_i), corrs.cells.co.shuffle.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,1);
        
        % b. solo vs. co norm max
        co_allo_cell_mat=co_allo_cell_mat./max(co_allo_cell_mat,[],2);
        solo_allo_cell_mat=solo_allo_cell_mat./max(solo_allo_cell_mat,2);
        [corrs.norm_PV.co.shuffle.mat(:,ego_bin_i,shuffle_i), corrs.norm_PV.co.shuffle.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
        
        %c. other dir solo vs. co
        co_allo_cell_mat=PSTH_solo_shuffle_ego_bins(:,:,shuffle_i,ego_bin_i);
        solo_allo_cell_mat=solo_PSTH_other_dir;
        [corrs.PV.other_dir.shuffle.mat(:,ego_bin_i,shuffle_i), corrs.PV.other_dir.shuffle.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
        [corrs.cells.other_dir.shuffle.mat(:,ego_bin_i,shuffle_i), corrs.cells.other_dir.shuffle.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,1);
        
        % Lower bound - co permuted (shuffle) vs. solo (all)
        %=====================================================
        % a. solo vs co
        co_allo_cell_mat_not_shuffled=PSTH_co_ego_bins(:,:,ego_bin_i);
        solo_allo_cell_mat=solo_PSTH;
        %co_allo_cell_mat=co_allo_cell_mat(randperm(size(co_allo_cell_mat_not_shuffled,1)),randperm(size(co_allo_cell_mat_not_shuffled,2)));
        
        %permute positions (rows)
        co_allo_cell_mat=co_allo_cell_mat_not_shuffled(randperm(size(co_allo_cell_mat_not_shuffled,1)),:);
        [corrs.PV.co.lower_bound.mat(:,ego_bin_i,shuffle_i), corrs.PV.co.lower_bound.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
        
        co_allo_cell_mat=co_allo_cell_mat./max(co_allo_cell_mat,[],2);
        solo_allo_cell_mat=solo_allo_cell_mat./max(solo_allo_cell_mat,2);
        [corrs.norm_PV.co.lower_bound.mat(:,ego_bin_i,shuffle_i), corrs.norm_PV.co.lower_bound.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
       
        %permute cells (columns)
        co_allo_cell_mat=co_allo_cell_mat_not_shuffled(:,randperm(size(co_allo_cell_mat_not_shuffled,2)));
        solo_allo_cell_mat=solo_PSTH;
        [corrs.cells.co.lower_bound.mat(:,ego_bin_i,shuffle_i), corrs.cells.co.lower_bound.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,1);
        
        % a. solo vs co
       
       
%         %c. other dir solo vs. co - norm to max
%         co_allo_cell_mat=co_allo_cell_mat./max(co_allo_cell_mat,[],2);
%         solo_allo_cell_mat=solo_allo_cell_mat./max(solo_allo_cell_mat,2);
%         [corrs.norm_PV.other_dir.shuffle.mat(:,ego_bin_i,shuffle_i), corrs.norm_PV.other_dir.shuffle.tuning(ego_bin_i,shuffle_i)]=corr_mat_cells_by_pos(co_allo_cell_mat,solo_allo_cell_mat,2);
%         
    end
    
end

