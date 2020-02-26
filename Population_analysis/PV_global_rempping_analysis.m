% population vector - global remapping analysis
clear
% parameters:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_info=dir(dir_data);
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
load(fullfile(param_folder,'solo_params.mat'));
pop_vec_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\pop_vec\';

%shuffle params:
n_shuffles=10;
%cell selection:
SI_threshold=1;
min_n_spike=50;
solo_time_spent_minimum_for_1D_bins=0.03;%need to check!
co_bins=-40:10:40;

%figure prop:
figure('units','normalized','outerposition',[0 0 1 1])
%pos:
width=0.15;
height=0.15;
hor_dis=0.5;
ver_dis=0.25;

x_pos(1)=0.1;
x_pos(2)=x_pos(1)+hor_dis;
%x_pos(3)=x_pos(2)+hor_dis;

y_pos(1)=0.8;
y_pos(2)=y_pos(1)-ver_dis;
y_pos(3)=y_pos(2)-ver_dis;

%% 1. load all data:
for dir_i=1:2
    cell_count=0;
    
    PV_co_dir_bin=[];PV_solo_dir=[];
    for cell_i=3:length({dir_info.name})-1
        % load data:
        load(fullfile(dir_data,dir_info(cell_i).name))
        num_spike_during_flight=sum(isnan(cell_co_solo_initial_analysis.solo(1).spikes.ts_usec(:)))+sum(isnan(cell_co_solo_initial_analysis.co(1).spikes.ts_usec(:)));
        if num_spike_during_flight<=min_n_spike || cell_co_solo_initial_analysis.solo(dir_i).SI<SI_threshold
            continue
        end
        cell_count=cell_count+1;
        
        %load flight solo data:
        solo_flight_bsp_ts=cell_co_solo_initial_analysis.solo(dir_i).bsp.ts_usec;
        solo_flight_spike_ts=cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec;
        solo_flight_pos_bsp=cell_co_solo_initial_analysis.solo(dir_i).bsp.x_pos;
        solo_flight_pos_spike=cell_co_solo_initial_analysis.solo(dir_i).spikes.x_pos;
        
        %load solo data:
        solo_pos_bsp=cell_co_solo_initial_analysis.solo(dir_i).bsp.x_pos';solo_pos_bsp=solo_pos_bsp(:);
        solo_pos_spike=cell_co_solo_initial_analysis.solo(dir_i).spikes.x_pos';solo_pos_spike=solo_pos_spike(:);
        %remove nans
        solo_pos_bsp=solo_pos_bsp(~isnan(solo_pos_bsp));
        solo_pos_spike=solo_pos_spike(~isnan(solo_pos_spike));
        solo_tuning_curve=cell_co_solo_initial_analysis.solo(dir_i).x_pos_firing_rate{1, 1};
        %load co data:
        co_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.x_pos';co_pos_bsp=co_pos_bsp(:);
        co_pos_spike=cell_co_solo_initial_analysis.co(dir_i).spikes.x_pos';co_pos_spike=co_pos_spike(:);
        co_dis_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.dis_m';co_dis_bsp=co_dis_bsp(:);
        co_dis_spike=cell_co_solo_initial_analysis.co(dir_i).spikes.dis_m';co_dis_spike=co_dis_spike(:);
        n_co=size(cell_co_solo_initial_analysis.co(dir_i).bsp.dis_m,1);
        %remove nans
        co_pos_bsp=co_pos_bsp(~isnan(co_pos_bsp));
        co_pos_spike=co_pos_spike(~isnan(co_pos_spike));
        co_dis_bsp=co_dis_bsp(~isnan(co_dis_bsp));
        co_dis_spike=co_dis_spike(~isnan(co_dis_spike));
        
        %         %% plot
        %             subplot(3,8,1)
        %             plot(solo_X_bins_vector_of_centers,solo_tuning_curve)
        %
        %             subplot(3,8,2)
        %             plot(solo_X_bins_vector_of_centers,cell_co_solo_initial_analysis.co(dir_i).firing_rate.allo_x_pos(1,:))
        %
        %% comupte tuning curves fo co data bin by bin of CO:
        for bin_i=1:length(co_bins)-1
            relevant_ind_bsp=(co_dis_bsp>=co_bins(bin_i) & co_dis_bsp<co_bins(bin_i+1));
            relevant_ind_spike=(co_dis_spike>=co_bins(bin_i) & co_dis_spike<co_bins(bin_i+1));
            relevant_shuffle_size_per_flight=round(sum(relevant_ind_bsp)/n_co);
            co_i_bsp_pos=co_pos_bsp(relevant_ind_bsp);
            co_i_spike_pos=co_pos_spike(relevant_ind_spike);
            if ~isempty(co_i_spike_pos)
                [~, ~, ~, PV_co_dir_bin(:,cell_count,bin_i), ~,~] ...
                    = fn_compute_generic_1D_tuning_new_smooth ...
                    (co_i_bsp_pos, co_i_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                %                   %% plot
                %             subplot(3,8,bin_i+8)
                %             plot(solo_X_bins_vector_of_centers,PV_co_dir_bin(:,cell_count,bin_i))
            end
            
            %% create shuffle
            for shuffle_i=1:n_shuffles
                bsp_data_flight_to_remove_ind_all=[];
                spike_data_flight_to_remove_ind_all=[];
                
                [ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(find(relevant_ind_bsp),length(relevant_ind_bsp));
                for epoch_i=1:length(ind_length)% run over all CO that had data in this distance bin
                    %find the epoch positions
                    epoch_pos=co_pos_bsp(start_ind(epoch_i):end_ind(epoch_i));
                    % find solo data flight with this positions
                    [row,col]=find(solo_flight_pos_bsp>=epoch_pos(1) & solo_flight_pos_bsp<=epoch_pos(end));
                    if isempty(row)
                        continue
                    end
                    relevanlt_flights=unique(row);
                    %choose random flight:
                    flight_to_remove_data_ind=randi(length(relevanlt_flights));
                    flight_chosen=relevanlt_flights(flight_to_remove_data_ind);
                    %take the data to shuffle i bsp:
                    len_epoch_to_remove=length(col(row==flight_chosen));
                    bsp_data_flight_to_remove_ind=[repmat(flight_chosen,len_epoch_to_remove,1),col(row==flight_chosen)];
                    %take the data to shuffle i spikes:
                    spikes_pos_flight_i=solo_flight_pos_spike(flight_chosen,:);
                    ind_spike_to_remove=find(spikes_pos_flight_i()>=epoch_pos(1) & spikes_pos_flight_i<=epoch_pos(end));
                    if ~isempty(ind_spike_to_remove)
                        spike_data_flight_to_remove_ind=[repmat(flight_chosen,length(ind_spike_to_remove),1),ind_spike_to_remove'];
                    else
                        spike_data_flight_to_remove_ind=[];
                    end
                    bsp_data_flight_to_remove_ind_all=[bsp_data_flight_to_remove_ind_all;bsp_data_flight_to_remove_ind];
                    spike_data_flight_to_remove_ind_all=[spike_data_flight_to_remove_ind_all;spike_data_flight_to_remove_ind];
                    
                end
                if ~isempty(spike_data_flight_to_remove_ind_all)
                    %find the spikes ans data of shuffle and all the
                    %rest:
                    remove_data_bsp_pos=solo_flight_pos_bsp(sub2ind(size(solo_flight_pos_bsp),bsp_data_flight_to_remove_ind_all(:,1),bsp_data_flight_to_remove_ind_all(:,2)));
                    remove_data_spike_pos=solo_flight_pos_spike(sub2ind(size(solo_flight_pos_spike),spike_data_flight_to_remove_ind_all(:,1),spike_data_flight_to_remove_ind_all(:,2)));
                    rest_of_data_bsp_pos=solo_flight_pos_bsp;
                    rest_of_data_bsp_pos(sub2ind(size(solo_flight_pos_bsp),bsp_data_flight_to_remove_ind_all(:,1),bsp_data_flight_to_remove_ind_all(:,2)))=nan;
                    rest_of_data_bsp_pos=rest_of_data_bsp_pos(:);
                    rest_of_data_spike_pos=solo_flight_pos_spike;
                    rest_of_data_spike_pos(sub2ind(size(solo_flight_pos_spike),spike_data_flight_to_remove_ind_all(:,1),spike_data_flight_to_remove_ind_all(:,2)))=nan;
                    rest_of_data_spike_pos=rest_of_data_spike_pos(:);
                    %remove nans:
                    remove_data_bsp_pos=remove_data_bsp_pos(~isnan(remove_data_bsp_pos));
                    remove_data_spike_pos=remove_data_spike_pos(~isnan(remove_data_spike_pos));
                    rest_of_data_bsp_pos=rest_of_data_bsp_pos(~isnan(rest_of_data_bsp_pos));
                    rest_of_data_spike_pos=rest_of_data_spike_pos(~isnan(rest_of_data_spike_pos));
                    % compute remove data tuning curve
                    [~, ~, ~, remove_data_tuning_curve(:,cell_count,shuffle_i,bin_i), ~,~] ...
                        = fn_compute_generic_1D_tuning_new_smooth ...
                        (remove_data_bsp_pos, remove_data_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                    
                    % compute rest of data tuning curve
                    [~, ~, ~, rest_of_data_tuning_curve(:,cell_count,shuffle_i,bin_i), ~,~] ...
                        = fn_compute_generic_1D_tuning_new_smooth ...
                        (rest_of_data_bsp_pos, rest_of_data_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                else
                    remove_data_tuning_curve(:,cell_count,shuffle_i,bin_i)=nan*zeros(length(solo_X_bins_vector_of_centers),1);
                    rest_of_data_tuning_curve(:,cell_count,shuffle_i,bin_i)=nan*zeros(length(solo_X_bins_vector_of_centers),1);
                end
            end
            %             n_epoch=min(n_co,size(solo_flight_bsp_ts,1));
            %             if n_epoch<=size(solo_flight_bsp_ts,1)
            %                 flight_vec=1:n_epoch;
            %             else
            %                 flight_vec=randi(size(solo_flight_bsp_ts,1),n_epoch,1);
            %             end
            %             remove_data_bsp_pos=[];
            %             remove_data_spike_pos=[];
            %             rest_of_data_bsp_pos=[];
            %             rest_of_data_spike_pos=[];
            %
            %             for flight_i=1:size(solo_flight_bsp_ts,1)
            %                 if ismember(flight_i,flight_vec)
            %                     %find non nan ind for flight
            %                     flight_non_nan_vec=1:max(find(~isnan(solo_flight_bsp_ts(flight_i,:))));
            %                     len_epoch_flight_to_remove=min(relevant_shuffle_size_per_flight,length(flight_non_nan_vec));
            %                    remove_bsp_ind=zeros(length(solo_flight_bsp_ts(flight_i,:)),1);
            %
            %                     if length(flight_non_nan_vec)<=relevant_shuffle_size_per_flight+1
            %                         remove_bsp_ind(1:len_epoch_flight_to_remove)=1;
            %                     else
            %                         remove_data_data_start=randi(length(flight_non_nan_vec)-len_epoch_flight_to_remove-1);
            %                         remove_bsp_ind(remove_data_data_start:remove_data_data_start+len_epoch_flight_to_remove)=1;
            %
            %                     end
            %                     remove_bsp_ind=logical(remove_bsp_ind);
            %                     remove_bsp_data_ts=solo_flight_bsp_ts(flight_i,remove_bsp_ind);
            %                     flight_spike_ts=solo_flight_spike_ts(flight_i,:);
            %                     remove_spike_ind=(flight_spike_ts>=remove_bsp_data_ts(1) & flight_spike_ts<=remove_bsp_data_ts(end));
            %
            %                     %
            %                     remove_data_bsp_pos=[remove_data_bsp_pos solo_flight_pos_bsp(flight_i,remove_bsp_ind)];
            %                     remove_data_spike_pos=[remove_data_spike_pos solo_flight_pos_spike(flight_i,remove_spike_ind)];
            %                     rest_of_data_bsp_pos=[rest_of_data_bsp_pos solo_flight_pos_bsp(flight_i,~remove_bsp_ind)];
            %                     rest_of_data_spike_pos=[rest_of_data_spike_pos solo_flight_pos_spike(flight_i,~remove_spike_ind)];
            %                 else
            %                     rest_of_data_bsp_pos=[rest_of_data_bsp_pos solo_flight_pos_bsp(flight_i,:)];
            %                     rest_of_data_spike_pos=[rest_of_data_spike_pos solo_flight_pos_spike(flight_i,:)];
            %                 end
            %
            %             end
            %             %remove nans
            %             rest_of_data_bsp_pos=rest_of_data_bsp_pos(~isnan(rest_of_data_bsp_pos));
            %             rest_of_data_spike_pos=rest_of_data_spike_pos(~isnan(rest_of_data_spike_pos));
            %             remove_data_bsp_pos=remove_data_bsp_pos(~isnan(remove_data_bsp_pos));
            %             remove_data_spike_pos=remove_data_spike_pos(~isnan(remove_data_spike_pos));
            %
            %             % compute remove data tuning curve
            %             [~, ~, ~, remove_data_tuning_curve(:,cell_count,shuffle_i,bin_i), ~,~] ...
            %                 = fn_compute_generic_1D_tuning_new_smooth ...
            %                 (remove_data_bsp_pos, remove_data_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
            %
            %             % compute rest of data tuning curve
            %             [~, ~, ~, rest_of_data_tuning_curve(:,cell_count,shuffle_i,bin_i), ~,~] ...
            %                 = fn_compute_generic_1D_tuning_new_smooth ...
            %                 (rest_of_data_bsp_pos, rest_of_data_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
            %
            %         end
            %          %% plot
            %             subplot(3,8,bin_i+16)
            %             plot(solo_X_bins_vector_of_centers,squeeze(rest_of_data_tuning_curve(:,cell_count,:,bin_i)))
        end
        PV_solo_dir(:,cell_count)=solo_tuning_curve;
        
        
        %      % save cell figure
        %      fig_name=fullfile(pop_vec_folder,['pop_vector_fig_cell_',num2str(cell_co_solo_initial_analysis.exp_data.cell_num),'_dir',num2str(dir_i),'.png']);
        %             saveas(gcf,fig_name)
        %             clf
    end
    norm_PV_solo_dir=PV_solo_dir./max(PV_solo_dir);
    
    for bin_i=1:length(co_bins)-1
        PV_co_dir_bin_i=PV_co_dir_bin(:,:,bin_i);
        norm_PV_co_dir_bin_i=PV_co_dir_bin_i./max(PV_co_dir_bin_i);
        PV_corr=corr(norm_PV_solo_dir',norm_PV_co_dir_bin_i','rows','pairwise');
        cell_corr=corr(norm_PV_solo_dir,norm_PV_co_dir_bin_i,'rows','pairwise');
        cell_corr_bin(:,bin_i)=diag(cell_corr);
        cell_corr_mean_bin(bin_i)=nanmean(cell_corr_bin(:,bin_i));
        PV_corr_bin(:,bin_i)=diag(PV_corr);
        corr_bin(bin_i)=nanmean(PV_corr_bin(:,bin_i));
        for shuffle_i=1:n_shuffles
            PV_solo_remove_shuffle=remove_data_tuning_curve(:,:,shuffle_i,bin_i);
            PV_solo_rest_of_data_shuffle=rest_of_data_tuning_curve(:,:,shuffle_i,bin_i);
            norm_PV_solo_remove_shuffle=PV_solo_remove_shuffle./max(PV_solo_remove_shuffle);
            norm_PV_solo_rest_of_data_shuffle=PV_solo_rest_of_data_shuffle./max(PV_solo_rest_of_data_shuffle);
            PV_corr_shuffle=corr(norm_PV_solo_rest_of_data_shuffle',norm_PV_solo_remove_shuffle','rows','pairwise');
            cell_corr_shuffle=corr(norm_PV_solo_rest_of_data_shuffle,norm_PV_solo_remove_shuffle,'rows','pairwise');
            cell_corr_shuffle_i(:,shuffle_i,bin_i)=diag(cell_corr_shuffle);
            cell_corr_mean_shuffle(shuffle_i,bin_i)=nanmean(cell_corr_shuffle_i(:,shuffle_i,bin_i));
            PV_corr_shuffle(:,shuffle_i,bin_i)=diag(PV_corr_shuffle);
            corr_shuffle(shuffle_i,bin_i)=nanmean(PV_corr_shuffle(:,shuffle_i,bin_i));
        end
    end
    
    %% plots
    %1. plot hist of cells corr
    x_vec=-1:0.1:1;
    axes('position',[x_pos(dir_i),y_pos(1),width,height])
    h=hist(nanmean(cell_corr_shuffle_i(:,4:5),2),x_vec);
    b1 = bar(x_vec,h/sum(h),'FaceColor',[0.5,0.5,0.5]);
    set(get(b1,'Children'),'FaceAlpha',0.3)
    hold on
    h=hist(nanmean(cell_corr_bin(:,4:5),2),x_vec);
    b2 = bar(x_vec,h/sum(h),'FaceColor',[255,0,255]./255);
    set(get(b2,'Children'),'FaceAlpha',0.2)
    title('Cell correlations')
    l=legend('shuffle (solo vs solo)','solo vs co (only in +-10)','location','eastoutside');
    set(l, 'position',[x_pos(dir_i)+1.05*width y_pos(1) 0.1 0.05])
    xlabel('r')
    ylabel('proportion of cells')
    % plot cell cor by bin vs shuffle
    axes('position',[x_pos(dir_i),y_pos(2),width,height])
    plot(co_bins(1:end-1)+5,cell_corr_mean_shuffle,'color',[.5 .5 .5])
    hold on;
    plot(co_bins(1:end-1)+5,cell_corr_mean_bin,'linewidth',2,'color','m')
    ylabel('r')
    xlabel('Inter-bat distance (m)')
    title('cell corr')
    % plot cell cor by bin vs shuffle
    axes('position',[x_pos(dir_i),y_pos(3),width,height])
    plot(co_bins(1:end-1)+5,corr_shuffle,'color',[.5 .5 .5])
    hold on;
    plot(co_bins(1:end-1)+5,corr_bin,'linewidth',2,'color','m')
    ylabel('r')
    xlabel('Inter-bat distance (m)')
    title('Population vector corr')
    
    
    
end
fig_name=fullfile(pop_vec_folder,['pop_vector_fig','.png']);
saveas(gcf,fig_name)