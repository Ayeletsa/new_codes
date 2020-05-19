% population vector - global remapping analysis
clear
% parameters:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_info=dir(dir_data);
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
load(fullfile(param_folder,'solo_params.mat'));
pop_vec_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\pop_vec\';


%p = parpool(4);
%shuffle params:
n_shuffles=100;
%cell selection:
SI_threshold=1;
min_n_spike_vec=[100];
solo_time_spent_minimum_for_1D_bins=0.03;%need to check!
range=[-40 40];% change to run just on non-overlap part between solo and co
window=10;
step=2;
bin_start=range(1):step:range(2)-window;
bin_end=range(1)+window:step:range(2);
bin_center=[bin_start+bin_end]./2;
bin_dis=range(1):step:range(2);
run_only_on_intersent_data=0;
full_data=0;
%solo_X_n_bins = solo_X_max / 2;
solo_X_bin_size_vec = [1];
for solo_X_bin_size=solo_X_bin_size_vec
    
    
    for min_n_spike=min_n_spike_vec
        %% 1. load all data:
        for dir_i=1:2
            clearvars -except  full_data solo_X_max solo_X_min solo_X_bin_size dir_i min_n_spike solo_X_bin_size_vec min_n_spike_vec run_only_on_intersent_data solo_X_bins_vector_of_centers solo_X_bins_vector solo_X_bin_size bin_center bin_end bin_start step window range solo_time_spent_minimum_for_1D_bins min_n_spike dir_data param_folder dir_info dir_i SI_threshold n_shuffles pop_vec_folder dir_data bin_dis
           
            bins_ii=solo_X_bin_size;
            load(fullfile(param_folder,'solo_params.mat'));
            solo_X_bin_size=bins_ii;
            solo_X_bins_vector=solo_X_min:solo_X_bin_size:solo_X_max;
            solo_X_bins_vector_of_centers=solo_X_bins_vector(1:end-1)+solo_X_bin_size/2;
            solo_time_spent_minimum_for_1D_bins=0.03;
            
            cell_count=0;
            
            PV_co_dir_bin=[];PV_solo_dir=[];
            day=0;
            day_i=0;
            file_names={dir_info.name};
            file_names=file_names(find([dir_info.isdir]==0));
            %%
            for cell_i=1:length(file_names)
                % for cell_i=3:15
                
                % load data:
                load(fullfile(dir_data,file_names{cell_i}))
                num_spike_during_flight=sum(~isnan(cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(dir_i).spikes.ts_usec(:)));
                if num_spike_during_flight<=min_n_spike || cell_co_solo_initial_analysis.solo(dir_i).SI<SI_threshold
                    continue
                end
                cell_count=cell_count+1;
                prev_day=day;
                day=cell_co_solo_initial_analysis.exp_data.day;
                
                %load flight solo data:
                solo_flight_bsp_ts=cell_co_solo_initial_analysis.solo(dir_i).bsp.ts_usec;
                solo_flight_x_pos_bsp=cell_co_solo_initial_analysis.solo(dir_i).bsp.x_pos;
                solo_flight_y_pos_bsp=cell_co_solo_initial_analysis.solo(dir_i).bsp.y_pos;
                solo_flight_x_pos_spike=cell_co_solo_initial_analysis.solo(dir_i).spikes.x_pos;
                solo_flight_y_pos_spike=cell_co_solo_initial_analysis.solo(dir_i).spikes.y_pos;
                solo_flight_spike_ts=cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec;
                
                %load co data:
                co_flight_bsp_ts=cell_co_solo_initial_analysis.co(dir_i).bsp.ts_usec;
                co_flight_x_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.x_pos;
                co_flight_y_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.y_pos;
                co_flight_x_pos_spike=cell_co_solo_initial_analysis.co(dir_i).spikes.x_pos;
                co_flight_y_pos_spike=cell_co_solo_initial_analysis.co(dir_i).spikes.y_pos;
                co_flight_dis_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.dis_m;
                co_flight_dis_spike=cell_co_solo_initial_analysis.co(dir_i).spikes.dis_m;
                
                %other dir:
                other_dir_solo_flight_x_pos_bsp=cell_co_solo_initial_analysis.solo(3-dir_i).bsp.x_pos;
                other_dir_solo_flight_x_pos_spike=cell_co_solo_initial_analysis.solo(3-dir_i).spikes.x_pos;
                other_dir_solo_flight_y_pos_bsp=cell_co_solo_initial_analysis.solo(3-dir_i).bsp.y_pos;
                other_dir_solo_flight_y_pos_spike=cell_co_solo_initial_analysis.solo(3-dir_i).spikes.y_pos;
                
                
                if run_only_on_intersent_data==1
                    %load flight solo data:
                    copy_solo_flight_bsp_ts=solo_flight_bsp_ts;
                    copy_solo_flight_x_pos_bsp=solo_flight_x_pos_bsp;
                    copy_solo_flight_y_pos_bsp=solo_flight_y_pos_bsp;
                    copy_solo_flight_x_pos_spike=solo_flight_x_pos_spike;
                    copy_solo_flight_y_pos_spike=solo_flight_y_pos_spike;
                    copy_solo_flight_spike_ts=solo_flight_spike_ts;
                    
                    %load co data:
                    copy_co_flight_bsp_ts=co_flight_bsp_ts;
                    copy_co_flight_x_pos_bsp=co_flight_x_pos_bsp;
                    copy_co_flight_y_pos_bsp=co_flight_y_pos_bsp;
                    copy_co_flight_x_pos_spike=co_flight_x_pos_spike;
                    copy_co_flight_y_pos_spike=co_flight_y_pos_spike;
                    copy_co_flight_dis_bsp=co_flight_dis_bsp;
                    copy_co_flight_dis_spike=co_flight_dis_spike;
                    
                    solo_flight_bsp_ts=nan*zeros(size(solo_flight_bsp_ts));
                    solo_flight_x_pos_bsp=nan*zeros(size(solo_flight_bsp_ts));
                    solo_flight_y_pos_bsp=nan*zeros(size(solo_flight_bsp_ts));
                    solo_flight_x_pos_spike=nan*zeros(size(solo_flight_x_pos_spike));
                    solo_flight_y_pos_spike=nan*zeros(size(solo_flight_x_pos_spike));
                    solo_flight_spike_ts=nan*zeros(size(solo_flight_x_pos_spike));
                    
                    %load co data:
                    co_flight_bsp_ts=nan*zeros(size(co_flight_x_pos_bsp));
                    co_flight_x_pos_bsp=nan*zeros(size(co_flight_x_pos_bsp));
                    co_flight_y_pos_bsp=nan*zeros(size(co_flight_x_pos_bsp));
                    co_flight_x_pos_spike=nan*zeros(size(co_flight_x_pos_spike));
                    co_flight_y_pos_spike=nan*zeros(size(co_flight_x_pos_spike));
                    co_flight_dis_bsp=nan*zeros(size(co_flight_x_pos_bsp));
                    co_flight_dis_spike=nan*zeros(size(co_flight_x_pos_spike));
                    
                    %other dir:
                    copy_other_dir_solo_flight_x_pos_bsp=other_dir_solo_flight_x_pos_bsp;
                    copy_other_dir_solo_flight_x_pos_spike=other_dir_solo_flight_x_pos_spike;
                    copy_other_dir_solo_flight_y_pos_bsp=other_dir_solo_flight_y_pos_bsp;
                    copy_other_dir_solo_flight_y_pos_spike=other_dir_solo_flight_y_pos_spike;
                    
                    
                    for bin_i=1:length(solo_X_bins_vector_of_centers)
                        solo_x_relevant_bsp_ind=(copy_solo_flight_x_pos_bsp>=solo_X_bins_vector(bin_i)& copy_solo_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                        if sum(solo_x_relevant_bsp_ind(:))~=0
                            solo_y_bsp_pos_per_x_bin=copy_solo_flight_y_pos_bsp(solo_x_relevant_bsp_ind);
                            solo_10_perc(bin_i)=prctile(solo_y_bsp_pos_per_x_bin,10);
                            solo_90_perc(bin_i)=prctile(solo_y_bsp_pos_per_x_bin,90);
                            
                            co_relevant_bsp_ind=(copy_co_flight_x_pos_bsp>=solo_X_bins_vector(bin_i)& copy_co_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                            co_y_bsp_pos_per_x_bin=copy_co_flight_y_pos_bsp(co_relevant_bsp_ind);
                            co_10_perc(bin_i)=prctile(co_y_bsp_pos_per_x_bin,10);
                            co_90_perc(bin_i)=prctile(co_y_bsp_pos_per_x_bin,90);
                            
                            if sum(isnan([co_10_perc(bin_i),solo_10_perc(bin_i),co_90_perc(bin_i),solo_90_perc(bin_i)]))==0
                                range_intersect=[max([co_10_perc(bin_i),solo_10_perc(bin_i)]), min([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                                range_union=[min([co_10_perc(bin_i),solo_10_perc(bin_i)]), max([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                                intersect_by_union_cell(bin_i)=abs(range_intersect(2)-range_intersect(1))/abs(range_union(2)-range_union(1));
                                %check if there is no intersection
                                if solo_10_perc(bin_i)>co_90_perc(bin_i) | co_10_perc(bin_i)>solo_90_perc(bin_i)
                                    intersect_by_union_cell(bin_i)=0;
                                    range_intersect=[nan nan];
                                end
                            else
                                range_intersect=[max([co_10_perc(bin_i),solo_10_perc(bin_i)]), min([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                                intersect_by_union_cell(bin_i)=nan;
                            end
                            if full_data==1
                                range_intersect=[nanmax([min(solo_y_bsp_pos_per_x_bin),min(co_y_bsp_pos_per_x_bin)]),nanmin([max(solo_y_bsp_pos_per_x_bin),max(co_y_bsp_pos_per_x_bin)])];
                                %check if there is no intersection
                                if min(solo_y_bsp_pos_per_x_bin)>max(co_y_bsp_pos_per_x_bin) | min(co_y_bsp_pos_per_x_bin)>max(solo_y_bsp_pos_per_x_bin)
                                    intersect_by_union_cell(bin_i)=0;
                                    range_intersect=[nan nan];
                                end
                            end
                            
                            
                            %remove data out of range:
                            %1. for solo bsp data:
                            solo_bsp_ind_for_analysis=zeros(size(copy_solo_flight_y_pos_bsp));
                            solo_bsp_y_relevant_ind=(copy_solo_flight_y_pos_bsp>=range_intersect(1) & copy_solo_flight_y_pos_bsp<=range_intersect(2));
                            solo_bsp_ind_for_analysis=(solo_x_relevant_bsp_ind==1 & solo_bsp_y_relevant_ind==1);
                            solo_flight_bsp_ts(solo_bsp_ind_for_analysis)=copy_solo_flight_bsp_ts(solo_bsp_ind_for_analysis);
                            solo_flight_x_pos_bsp(solo_bsp_ind_for_analysis)=copy_solo_flight_x_pos_bsp(solo_bsp_ind_for_analysis);
                            solo_flight_y_pos_bsp(solo_bsp_ind_for_analysis)=copy_solo_flight_y_pos_bsp(solo_bsp_ind_for_analysis);
                            
                            %2. for solo spike data:
                            solo_x_relevant_spike_ind=(copy_solo_flight_x_pos_spike>=solo_X_bins_vector(bin_i)& copy_solo_flight_x_pos_spike<solo_X_bins_vector(bin_i+1));
                            solo_spike_ind_for_analysis=zeros(size(solo_flight_y_pos_spike));
                            solo_spike_y_relevant_ind=(copy_solo_flight_y_pos_spike>=range_intersect(1) & copy_solo_flight_y_pos_spike<=range_intersect(2));
                            solo_spike_ind_for_analysis=(solo_x_relevant_spike_ind==1 & solo_spike_y_relevant_ind==1);
                            
                            solo_flight_spike_ts(solo_spike_ind_for_analysis)=copy_solo_flight_spike_ts(solo_spike_ind_for_analysis);
                            solo_flight_x_pos_spike(solo_spike_ind_for_analysis)=copy_solo_flight_x_pos_spike(solo_spike_ind_for_analysis);
                            solo_flight_y_pos_spike(solo_spike_ind_for_analysis)=copy_solo_flight_y_pos_spike(solo_spike_ind_for_analysis);
                            %3. for co bsp data:
                            co_x_relevant_bsp_ind=(copy_co_flight_x_pos_bsp>=solo_X_bins_vector(bin_i)& copy_co_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                            co_bsp_ind_for_analysis=zeros(size(copy_co_flight_y_pos_bsp));
                            co_bsp_y_relevant_ind=(copy_co_flight_y_pos_bsp>=range_intersect(1) & copy_co_flight_y_pos_bsp<=range_intersect(2));
                            co_bsp_ind_for_analysis=(co_x_relevant_bsp_ind==1 & co_bsp_y_relevant_ind==1);
                            
                            co_flight_x_pos_bsp(co_bsp_ind_for_analysis)=copy_co_flight_x_pos_bsp(co_bsp_ind_for_analysis);
                            co_flight_y_pos_bsp(co_bsp_ind_for_analysis)=copy_co_flight_y_pos_bsp(co_bsp_ind_for_analysis);
                            co_flight_dis_bsp(co_bsp_ind_for_analysis)=copy_co_flight_dis_bsp(co_bsp_ind_for_analysis);
                            co_flight_bsp_ts(co_bsp_ind_for_analysis)=copy_co_flight_bsp_ts(co_bsp_ind_for_analysis);
                            %4. for co spike data:
                            co_x_relevant_spike_ind=(copy_co_flight_x_pos_spike>=solo_X_bins_vector(bin_i)& copy_co_flight_x_pos_spike<solo_X_bins_vector(bin_i+1));
                            co_spike_ind_for_analysis=zeros(size(copy_co_flight_y_pos_spike));
                            co_spike_y_relevant_ind=(copy_co_flight_y_pos_spike>=range_intersect(1) & copy_co_flight_y_pos_spike<=range_intersect(2));
                            co_spike_ind_for_analysis=(co_x_relevant_spike_ind==1 & co_spike_y_relevant_ind==1);
                            
                            co_flight_x_pos_spike(co_spike_ind_for_analysis)=copy_co_flight_x_pos_spike(co_spike_ind_for_analysis);
                            co_flight_y_pos_spike(co_spike_ind_for_analysis)=copy_co_flight_y_pos_spike(co_spike_ind_for_analysis);
                            co_flight_dis_spike(co_spike_ind_for_analysis)=copy_co_flight_dis_spike(co_spike_ind_for_analysis);
                            
                            %5. for other dir bsp:
                            other_dir_x_relevant_bsp_ind=(copy_other_dir_solo_flight_x_pos_bsp>=solo_X_bins_vector(bin_i)& copy_other_dir_solo_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                            other_dir_bsp_ind_for_analysis=zeros(size(copy_co_flight_y_pos_spike));
                            other_dir_bsp_y_relevant_ind=(copy_other_dir_solo_flight_y_pos_bsp>=range_intersect(1) & copy_other_dir_solo_flight_y_pos_bsp<=range_intersect(2));
                            other_dir_bsp_ind_for_analysis=(other_dir_x_relevant_bsp_ind==1 & other_dir_bsp_y_relevant_ind==1);
                            
                            other_dir_solo_flight_x_pos_bsp(other_dir_bsp_ind_for_analysis)=copy_other_dir_solo_flight_x_pos_bsp(other_dir_bsp_ind_for_analysis);
                            other_dir_solo_flight_y_pos_bsp(other_dir_bsp_ind_for_analysis)=copy_other_dir_solo_flight_y_pos_bsp(other_dir_bsp_ind_for_analysis);
                            
                            %6. for other dir spike:
                            other_dir_x_relevant_spike_ind=(copy_other_dir_solo_flight_x_pos_spike>=solo_X_bins_vector(bin_i)& copy_other_dir_solo_flight_x_pos_spike<solo_X_bins_vector(bin_i+1));
                            other_dir_spike_ind_for_analysis=zeros(size(copy_co_flight_y_pos_spike));
                            other_dir_spike_y_relevant_ind=(copy_other_dir_solo_flight_y_pos_spike>=range_intersect(1) & copy_other_dir_solo_flight_y_pos_spike<=range_intersect(2));
                            other_dir_spike_ind_for_analysis=(other_dir_x_relevant_spike_ind==1 & other_dir_spike_y_relevant_ind==1);
                            
                            other_dir_solo_flight_x_pos_spike(other_dir_spike_ind_for_analysis)=copy_other_dir_solo_flight_x_pos_spike(other_dir_spike_ind_for_analysis);
                            other_dir_solo_flight_y_pos_spike(other_dir_spike_ind_for_analysis)=copy_other_dir_solo_flight_y_pos_spike(other_dir_spike_ind_for_analysis);
                            
                        end
                    end
                end
                
                
                
                %load solo data:
                solo_x_pos_bsp=solo_flight_x_pos_bsp';solo_x_pos_bsp=solo_x_pos_bsp(:);
                solo_y_pos_bsp=solo_flight_y_pos_bsp';solo_y_pos_bsp=solo_y_pos_bsp(:);
                solo_x_pos_spike=solo_flight_x_pos_spike';solo_x_pos_spike=solo_x_pos_spike(:);
                solo_y_pos_spike=solo_flight_x_pos_spike';solo_y_pos_spike=solo_y_pos_spike(:);
                
                %remove nans
                solo_x_pos_bsp=solo_x_pos_bsp(~isnan(solo_x_pos_bsp));
                solo_y_pos_bsp=solo_y_pos_bsp(~isnan(solo_y_pos_bsp));
                solo_x_pos_spike=solo_x_pos_spike(~isnan(solo_x_pos_spike));
                solo_y_pos_spike=solo_y_pos_spike(~isnan(solo_y_pos_spike));
                %other dir:
                other_dir_solo_pos_bsp=other_dir_solo_flight_x_pos_bsp';other_dir_solo_pos_bsp=other_dir_solo_pos_bsp(:);
                other_dir_solo_pos_spike=other_dir_solo_flight_x_pos_spike';other_dir_solo_pos_spike=other_dir_solo_pos_spike(:);
                %remove nans
                other_dir_solo_pos_bsp=other_dir_solo_pos_bsp(~isnan(other_dir_solo_pos_bsp));
                other_dir_solo_pos_spike=other_dir_solo_pos_spike(~isnan(other_dir_solo_pos_spike));
                
                %load co data:
                co_x_pos_bsp=co_flight_x_pos_bsp';co_x_pos_bsp=co_x_pos_bsp(:);
                co_y_pos_bsp=co_flight_y_pos_bsp';co_y_pos_bsp=co_y_pos_bsp(:);
                co_y_pos_spike=co_flight_y_pos_spike';co_y_pos_spike=co_y_pos_spike(:);
                co_x_pos_spike=co_flight_x_pos_spike';co_x_pos_spike=co_x_pos_spike(:);
                co_dis_bsp=co_flight_dis_bsp';co_dis_bsp=co_dis_bsp(:);
                co_dis_spike=co_flight_dis_spike';co_dis_spike=co_dis_spike(:);
                n_co=size(cell_co_solo_initial_analysis.co(dir_i).bsp.dis_m,1);
                
                
                
                %remove nans
                co_x_pos_bsp=co_x_pos_bsp(~isnan(co_x_pos_bsp));
                co_y_pos_bsp=co_y_pos_bsp(~isnan(co_y_pos_bsp));
                co_x_pos_spike=co_x_pos_spike(~isnan(co_x_pos_spike));
                co_dis_bsp=co_dis_bsp(~isnan(co_dis_bsp));
                co_dis_spike=co_dis_spike(~isnan(co_dis_spike));
                co_y_pos_spike=co_y_pos_spike(~isnan(co_y_pos_spike));
                
                %         %% plot
                %             subplot(3,8,1)
                %             plot(solo_X_bins_vector_of_centers,solo_tuning_curve)
                %
                %             subplot(3,8,2)
                %             plot(solo_X_bins_vector_of_centers,cell_co_solo_initial_analysis.co(dir_i).firing_rate.allo_x_pos(1,:))
                if day~=prev_day
                    day_i=day_i+1;
                    % compute velocity:
                    %co_flight_y_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.y_pos;
                    %co_flight_x_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.x_pos;
                    %co_flight_bsp_ts=cell_co_solo_initial_analysis.co(dir_i).bsp.ts_usec;
                    %co_flight_dis=cell_co_solo_initial_analysis.co(dir_i).bsp.dis_m;
                    
                    idx=[];
                    for i=1:size(co_flight_dis_bsp,1)
                        [h ,a,b,idx(i,:)]=histcn(co_flight_dis_bsp(i,:)',bin_dis);
                        idx(i,find(idx(i,:)==0))=nan;
                    end
                    co_velocity=((sqrt(diff(co_flight_x_pos_bsp,[],2).^2+diff(co_flight_y_pos_bsp,[],2).^2)./diff(co_flight_bsp_ts,[],2)))*1e6;
                    co_velocity_plus_one=nan*zeros(size(co_velocity,1),size(co_velocity,2)+1);
                    co_velocity_plus_one(:,2:end)=co_velocity;
                    co_velocity_plus_one(:,1)=co_velocity_plus_one(:,2);
                    % compute mean y pos per x bin:
                    [h ,a,b,idx_x]=histcn(solo_x_pos_bsp,solo_X_bins_vector_of_centers);
                    co_flight_y_pos_bsp_x_bin=nan*zeros(size(co_flight_dis_bsp));
                    for bin_i=1:length(solo_X_bins_vector_of_centers)
                        mean_y_pos_per_bin=nanmean(solo_y_pos_bsp(find(idx_x==bin_i)));
                        relevant_ind=find(co_flight_x_pos_bsp>solo_X_bins_vector(bin_i) & co_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                        co_flight_y_pos_bsp_x_bin(relevant_ind)=abs(co_flight_y_pos_bsp(relevant_ind)-mean_y_pos_per_bin);
                        
                        solo_x_relevant_bsp_ind=(solo_flight_x_pos_bsp>=solo_X_bins_vector(bin_i)& solo_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                        if sum(solo_x_relevant_bsp_ind(:))~=0
                            solo_y_bsp_pos_per_x_bin=solo_flight_y_pos_bsp(solo_x_relevant_bsp_ind);
                            solo_10_perc(bin_i)=prctile(solo_y_bsp_pos_per_x_bin,10);
                            solo_90_perc(bin_i)=prctile(solo_y_bsp_pos_per_x_bin,90);
                            
                            co_relevant_bsp_ind=(co_flight_x_pos_bsp>=solo_X_bins_vector(bin_i)& co_flight_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                            co_y_bsp_pos_per_x_bin=co_flight_y_pos_bsp(co_relevant_bsp_ind);
                            co_10_perc(bin_i)=prctile(co_y_bsp_pos_per_x_bin,10);
                            co_90_perc(bin_i)=prctile(co_y_bsp_pos_per_x_bin,90);
                            
                            if sum(isnan([co_10_perc(bin_i),solo_10_perc(bin_i),co_90_perc(bin_i),solo_90_perc(bin_i)]))==0
                                range_intersect=[max([co_10_perc(bin_i),solo_10_perc(bin_i)]), min([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                                range_union=[min([co_10_perc(bin_i),solo_10_perc(bin_i)]), max([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                                intersect_by_union_cell(bin_i)=abs(range_intersect(2)-range_intersect(1))/abs(range_union(2)-range_union(1));
                                %fulldata
                                full_range_intersect=[nanmax([min(solo_y_bsp_pos_per_x_bin),min(co_y_bsp_pos_per_x_bin)]),nanmin([max(solo_y_bsp_pos_per_x_bin),max(co_y_bsp_pos_per_x_bin)])];
                                full_range_union=[nanmin([min(solo_y_bsp_pos_per_x_bin),min(co_y_bsp_pos_per_x_bin)]),nanmax([max(solo_y_bsp_pos_per_x_bin),max(co_y_bsp_pos_per_x_bin)])];
                                full_intersect_by_union_cell(bin_i)=abs(full_range_intersect(2)-full_range_intersect(1))/abs(full_range_union(2)-full_range_union(1));
                            else
                                range_intersect=[max([co_10_perc(bin_i),solo_10_perc(bin_i)]), min([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                                intersect_by_union_cell(bin_i)=nan;
                            end
                            if full_data==1
                                range_intersect=[nanmax([min(solo_y_bsp_pos_per_x_bin),min(co_y_bsp_pos_per_x_bin)]),nanmin([max(solo_y_bsp_pos_per_x_bin),max(co_y_bsp_pos_per_x_bin)])];
                            end
                        end
                    end
                    intersect_by_union_cell(intersect_by_union_cell==0)=[];
                    full_intersect_by_union_cell(full_intersect_by_union_cell==0)=[];
                    intersect_by_union(day_i)=nanmean(intersect_by_union_cell);
                    full_intersect_by_union(day_i)=nanmean(full_intersect_by_union_cell);
                    for ii=1:length(bin_dis)
                        [co, ind]=find(idx==ii);
                        vel_by_dis(day_i,ii)=nanmean(co_velocity_plus_one(sub2ind(size(co_velocity_plus_one),co, ind)));
                        y_pos_by_dis(day_i,ii)=nanmean(co_flight_y_pos_bsp_x_bin(sub2ind(size(co_velocity_plus_one),co, ind)));
                    end
                end
                %% comupte tuning curves for co data bin by bin of CO:
                for bin_i=1:length(bin_start)
                    bin_i
                    relevant_ind_bsp=(co_dis_bsp>=bin_start(bin_i) & co_dis_bsp<bin_end(bin_i));
                    relevant_ind_spike=(co_dis_spike>=bin_start(bin_i) & co_dis_spike<bin_end(bin_i));
                    relevant_shuffle_size_per_flight=round(sum(relevant_ind_bsp)/n_co);
                    co_i_bsp_pos=co_x_pos_bsp(relevant_ind_bsp);
                    co_i_spike_pos=co_x_pos_spike(relevant_ind_spike);
                    % if ~isempty(co_i_spike_pos)
                    [timespent_binned, ~, ~, PV_co_dir_bin(:,cell_count,bin_i), ~,~] ...
                        = fn_compute_generic_1D_tuning_new_smooth ...
                        (co_i_bsp_pos, co_i_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                    
                    %compute spatial information:
                    information_per_spike_co(cell_count,bin_i) = fn_compute_spatial_info (timespent_binned,PV_co_dir_bin(:,cell_count,bin_i));
                    
                    
                    %                   %% plot
                    %             subplot(3,8,bin_i+8)
                    %             plot(solo_X_bins_vector_of_centers,PV_co_dir_bin(:,cell_count,bin_i))
                    %   else
                    %  PV_co_dir_bin(:,cell_count,bin_i)=zeros(size(solo_X_bins_vector_of_centers));
                    % information_per_spike_co(cell_count,bin_i) = 0;
                    % end
                    [~, ~, ~, solo_tuning_curve, ~,~] ...
                        = fn_compute_generic_1D_tuning_new_smooth ...
                        (solo_x_pos_bsp, solo_x_pos_spike, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                    
                    %other dir tuning curve:
                    
                    [~, ~, ~, other_dir_solo_tuning_curve, ~,~] ...
                        = fn_compute_generic_1D_tuning_new_smooth ...
                        (other_dir_solo_pos_bsp, other_dir_solo_pos_spike, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                    
                    %% create shuffle uper bound (solo vs solo)
                    parfor shuffle_i=1:n_shuffles
                        shuffle_i
                        bsp_data_flight_to_remove_ind_all=[];
                        spike_data_flight_to_remove_ind_all=[];
                        
                        [ind_length,start_ind,end_ind]=find_length_of_consecutive_ind(find(relevant_ind_bsp),length(relevant_ind_bsp));
                        
                        for epoch_i=1:length(ind_length)% run over all CO that had data in this distance bin
                            
                            %find the epoch positions
                            epoch_pos=co_x_pos_bsp(start_ind(epoch_i):end_ind(epoch_i));
                            % find solo data flight with this positions
                            if dir_i==1
                                [row,col]=find(solo_flight_x_pos_bsp>=epoch_pos(1) & solo_flight_x_pos_bsp<=epoch_pos(end));
                            else
                                [row,col]=find(solo_flight_x_pos_bsp>=epoch_pos(end) & solo_flight_x_pos_bsp<=epoch_pos(1));
                                
                            end
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
                            spikes_pos_flight_i=solo_flight_x_pos_spike(flight_chosen,:);
                            if dir_i==1
                                
                                ind_spike_to_remove=find(spikes_pos_flight_i()>=epoch_pos(1) & spikes_pos_flight_i<=epoch_pos(end));
                            else
                                ind_spike_to_remove=find(spikes_pos_flight_i()>=epoch_pos(end) & spikes_pos_flight_i<=epoch_pos(1));
                                
                            end
                            if ~isempty(ind_spike_to_remove)
                                spike_data_flight_to_remove_ind=[repmat(flight_chosen,length(ind_spike_to_remove),1),ind_spike_to_remove'];
                            else
                                spike_data_flight_to_remove_ind=[];
                            end
                            bsp_data_flight_to_remove_ind_all=[bsp_data_flight_to_remove_ind_all;bsp_data_flight_to_remove_ind];
                            spike_data_flight_to_remove_ind_all=[spike_data_flight_to_remove_ind_all;spike_data_flight_to_remove_ind];
                            
                        end
                        %    if ~isempty(spike_data_flight_to_remove_ind_all)
                        %find the spikes ans data of shuffle and all the
                        %rest:
                        remove_data_bsp_pos=solo_flight_x_pos_bsp(sub2ind(size(solo_flight_x_pos_bsp),bsp_data_flight_to_remove_ind_all(:,1),bsp_data_flight_to_remove_ind_all(:,2)));
                        if ~isempty(spike_data_flight_to_remove_ind_all)
                            remove_data_spike_pos=solo_flight_x_pos_spike(sub2ind(size(solo_flight_x_pos_spike),spike_data_flight_to_remove_ind_all(:,1),spike_data_flight_to_remove_ind_all(:,2)));
                        else
                            remove_data_spike_pos=[];
                            
                        end
                        %rest_of_data_bsp_pos=solo_flight_pos_bsp;
                        %rest_of_data_bsp_pos(sub2ind(size(solo_flight_pos_bsp),bsp_data_flight_to_remove_ind_all(:,1),bsp_data_flight_to_remove_ind_all(:,2)))=nan;
                        %rest_of_data_bsp_pos=rest_of_data_bsp_pos(:);
                        %rest_of_data_spike_pos=solo_flight_pos_spike;
                        %if ~isemtpy()
                        %rest_of_data_spike_pos(sub2ind(size(solo_flight_pos_spike),spike_data_flight_to_remove_ind_all(:,1),spike_data_flight_to_remove_ind_all(:,2)))=nan;
                        %rest_of_data_spike_pos=rest_of_data_spike_pos(:);
                        %remove nans:
                        remove_data_bsp_pos=remove_data_bsp_pos(~isnan(remove_data_bsp_pos));
                        remove_data_spike_pos=remove_data_spike_pos(~isnan(remove_data_spike_pos));
                        %rest_of_data_bsp_pos=rest_of_data_bsp_pos(~isnan(rest_of_data_bsp_pos));
                        %rest_of_data_spike_pos=rest_of_data_spike_pos(~isnan(rest_of_data_spike_pos));
                        % compute remove data tuning curve
                        [timespent_binned, ~, ~, remove_data_tuning_curve(:,cell_count,shuffle_i,bin_i), ~,~] ...
                            = fn_compute_generic_1D_tuning_new_smooth ...
                            (remove_data_bsp_pos, remove_data_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                        %compute spatial information:
                        information_per_spike_shuffle(cell_count,shuffle_i,bin_i) = fn_compute_spatial_info (timespent_binned,remove_data_tuning_curve(:,cell_count,shuffle_i,bin_i));
                        
                        % compute rest of data tuning curve
                        %[~, ~, ~, rest_of_data_tuning_curve(:,cell_count,shuffle_i,bin_i), ~,~] ...
                        %   = fn_compute_generic_1D_tuning_new_smooth ...
                        %  (rest_of_data_bsp_pos, rest_of_data_spike_pos, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
                        % else
                        %    remove_data_tuning_curve(:,cell_count,shuffle_i,bin_i)=zeros(length(solo_X_bins_vector_of_centers),1);
                        %    information_per_spike_shuffle(cell_count,shuffle_i,bin_i) = 0;
                        %  rest_of_data_tuning_curve(:,cell_count,shuffle_i,bin_i)=zeros(length(solo_X_bins_vector_of_centers),1);
                        % end
                    end
                    
                    %          %% plot
                    %             subplot(3,8,bin_i+16)
                    %             plot(solo_X_bins_vector_of_centers,squeeze(rest_of_data_tuning_curve(:,cell_count,:,bin_i)))
                end
                PV_solo_dir(:,cell_count)=solo_tuning_curve;
                other_dir_PV_solo_dir(:,cell_count)=other_dir_solo_tuning_curve;
                
                %      % save cell figure
                %      fig_name=fullfile(pop_vec_folder,['pop_vector_fig_cell_',num2str(cell_co_solo_initial_analysis.exp_data.cell_num),'_dir',num2str(dir_i),'.png']);
                %             saveas(gcf,fig_name)
                %             clf
            end
            
            for bin_i=1:length(bin_start)
                bin_i
                
                PV_co_dir_bin_i=PV_co_dir_bin(:,:,bin_i);
                % normal
                [cell_corr_bin(:,bin_i) cell_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(PV_co_dir_bin_i,PV_solo_dir,1);
                [PV_corr_bin(:,bin_i) PV_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(PV_co_dir_bin_i,PV_solo_dir,2);
                
                %norm to max
                normMAX_PV_solo_dir=PV_solo_dir./(max(PV_solo_dir)+eps);
                normMAX_PV_co_dir_bin_i=PV_co_dir_bin_i./(max(PV_co_dir_bin_i)+eps);
                [normMAX_cell_corr_bin(:,bin_i) normMAX_cell_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_PV_solo_dir,1);
                [normMAX_PV_corr_bin(:,bin_i) normMAX_PV_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_PV_solo_dir,2);
                
                
                %norm to mean
                normMEAN_PV_solo_dir=PV_solo_dir./(nanmean(PV_solo_dir)+eps);
                normMEAN_PV_co_dir_bin_i=PV_co_dir_bin_i./(nanmean(PV_co_dir_bin_i)+eps);
                [normMEAN_cell_corr_bin(:,bin_i) normMEAN_cell_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_PV_solo_dir,1);
                [normMEAN_PV_corr_bin(:,bin_i) normMEAN_PV_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_PV_solo_dir,2);
                
                %lower bound of other dir:
                %================================================================
                % normal
                [other_dir_cell_corr_bin(:,bin_i) other_dir_cell_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(PV_co_dir_bin_i,other_dir_PV_solo_dir,1);
                [other_dir_PV_corr_bin(:,bin_i) other_dir_PV_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(PV_co_dir_bin_i,other_dir_PV_solo_dir,2);
                
                %norm to max
                normMAX_other_dir_PV_solo_dir=other_dir_PV_solo_dir./(max(other_dir_PV_solo_dir)+eps);
                normMAX_PV_co_dir_bin_i=PV_co_dir_bin_i./(max(PV_co_dir_bin_i)+eps);
                [other_dir_normMAX_cell_corr_bin(:,bin_i) other_dir_normMAX_cell_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_other_dir_PV_solo_dir,1);
                [other_dir_normMAX_PV_corr_bin(:,bin_i) other_dir_normMAX_PV_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_other_dir_PV_solo_dir,2);
                
                
                %norm to mean
                normMEAN_other_dir_PV_solo_dir=other_dir_PV_solo_dir./(nanmean(other_dir_PV_solo_dir)+eps);
                normMEAN_PV_co_dir_bin_i=PV_co_dir_bin_i./(nanmean(PV_co_dir_bin_i)+eps);
                [other_dir_normMEAN_cell_corr_bin(:,bin_i) other_dir_normMEAN_cell_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_other_dir_PV_solo_dir,1);
                [other_dir_normMEAN_PV_corr_bin(:,bin_i) other_dir_normMEAN_PV_corr_mean_bin(bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_other_dir_PV_solo_dir,2);
                
                
                parfor shuffle_i=1:n_shuffles
                    shuffle_i
                    % for upper bound shuffle (solo vs solo)
                    %=======================================================
                    PV_solo_remove_shuffle=remove_data_tuning_curve(:,:,shuffle_i,bin_i);
                    %PV_solo_rest_of_data_shuffle=rest_of_data_tuning_curve(:,:,shuffle_i,bin_i);
                    % normal
                    [cell_corr_bin_shuffle(:,shuffle_i,bin_i) cell_corr_mean_bin_shuffle(shuffle_i,bin_i)]=corr_mat_cells_by_pos(PV_solo_remove_shuffle,PV_solo_dir,1);
                    [PV_corr_bin_shuffle(:,shuffle_i,bin_i) PV_corr_mean_bin_shuffle(shuffle_i,bin_i)]=corr_mat_cells_by_pos(PV_solo_remove_shuffle,PV_solo_dir,2);
                    
                    %norm to max
                    normMAX_PV_solo_dir=PV_solo_dir./(max(PV_solo_dir)+eps);
                    normMAX_PV_co_dir_bin_i=PV_solo_remove_shuffle./(max(PV_solo_remove_shuffle)+eps);
                    [normMAX_cell_corr_bin_shuffle(:,shuffle_i,bin_i) normMAX_cell_corr_mean_bin_shuffle(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_PV_solo_dir,1);
                    [normMAX_PV_corr_bin_shuffle(:,shuffle_i,bin_i) normMAX_PV_corr_mean_bin_shuffle(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_PV_solo_dir,2);
                    
                    
                    %norm to mean
                    normMEAN_PV_solo_dir=PV_solo_dir./(nanmean(PV_solo_dir)+eps);
                    normMEAN_PV_co_dir_bin_i=PV_solo_remove_shuffle./(nanmean(PV_solo_remove_shuffle)+eps);
                    [normMEAN_cell_corr_bin_shuffle(:,shuffle_i,bin_i) normMEAN_cell_corr_mean_bin_shuffle(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_PV_solo_dir,1);
                    [normMEAN_PV_corr_bin_shuffle(:,shuffle_i,bin_i) normMEAN_PV_corr_mean_bin_shuffle(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_PV_solo_dir,2);
                    
                    % for lower bound shuffle (shuffle matrix)
                    %=======================================================
                    % shuffle rows:
                    PV_co_dir_bin_i_shuffle=PV_co_dir_bin_i(randperm(size(PV_co_dir_bin_i,1)),randperm(size(PV_co_dir_bin_i,2)));
                    % normal
                    [cell_corr_bin_shuffle_lower(:,shuffle_i,bin_i) cell_corr_mean_bin_shuffle_lower(shuffle_i,bin_i)]=corr_mat_cells_by_pos(PV_co_dir_bin_i_shuffle,PV_solo_dir,1);
                    [PV_corr_bin_shuffle_lower(:,shuffle_i,bin_i) PV_corr_mean_bin_shuffle_lower(shuffle_i,bin_i)]=corr_mat_cells_by_pos(PV_co_dir_bin_i_shuffle,PV_solo_dir,2);
                    
                    %norm to max
                    normMAX_PV_solo_dir=PV_solo_dir./(max(PV_solo_dir)+eps);
                    normMAX_PV_co_dir_bin_i=PV_co_dir_bin_i_shuffle./(max(PV_co_dir_bin_i_shuffle)+eps);
                    [normMAX_cell_corr_bin_shuffle_lower(:,shuffle_i,bin_i) normMAX_cell_corr_mean_bin_shuffle_lower(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_PV_solo_dir,1);
                    [normMAX_PV_corr_bin_shuffle_lower(:,shuffle_i,bin_i) normMAX_PV_corr_mean_bin_shuffle_lower(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMAX_PV_co_dir_bin_i,normMAX_PV_solo_dir,2);
                    
                    
                    %norm to mean
                    normMEAN_PV_solo_dir=PV_solo_dir./(nanmean(PV_solo_dir)+eps);
                    normMEAN_PV_co_dir_bin_i=PV_co_dir_bin_i_shuffle./(nanmean(PV_co_dir_bin_i_shuffle)+eps);
                    [normMEAN_cell_corr_bin_shuffle_lower(:,shuffle_i,bin_i) normMEAN_cell_corr_mean_bin_shuffle_lower(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_PV_solo_dir,1);
                    [normMEAN_PV_corr_bin_shuffle_lower(:,shuffle_i,bin_i) normMEAN_PV_corr_mean_bin_shuffle_lower(shuffle_i,bin_i)]=corr_mat_cells_by_pos(normMEAN_PV_co_dir_bin_i,normMEAN_PV_solo_dir,2);
                    
                    %             norm_PV_solo_remove_shuffle=PV_solo_remove_shuffle./(max(PV_solo_remove_shuffle));
                    %             norm_PV_solo_rest_of_data_shuffle=PV_solo_rest_of_data_shuffle./(max(PV_solo_rest_of_data_shuffle));
                    %             PV_corr_shuffle=corr(norm_PV_solo_rest_of_data_shuffle',norm_PV_solo_remove_shuffle','rows','pairwise');
                    %             cell_corr_shuffle=corr(norm_PV_solo_rest_of_data_shuffle,norm_PV_solo_remove_shuffle,'rows','pairwise');
                    %             cell_corr_shuffle_i(:,shuffle_i,bin_i)=diag(cell_corr_shuffle);
                    %             cell_corr_mean_shuffle(shuffle_i,bin_i)=nanmean(cell_corr_shuffle_i(:,shuffle_i,bin_i));
                    %             PV_corr_shuffle(:,shuffle_i,bin_i)=diag(PV_corr_shuffle);
                    %             corr_shuffle(shuffle_i,bin_i)=nanmean(PV_corr_shuffle(:,shuffle_i,bin_i));
                    %             %shuffle vs full solo data:
                    %             PV_corr_shuffle_full_solo=corr(norm_PV_solo_dir',norm_PV_solo_remove_shuffle','rows','pairwise');
                    %             cell_corr_shuffle_full_solo=corr(norm_PV_solo_dir,norm_PV_solo_remove_shuffle,'rows','pairwise');
                    %             cell_corr_shuffle_i_full_solo(:,shuffle_i,bin_i)=diag(cell_corr_shuffle_full_solo);
                    %             cell_corr_mean_shuffle_full_solo(shuffle_i,bin_i)=nanmean(cell_corr_shuffle_i_full_solo(:,shuffle_i,bin_i));
                    %             PV_corr_shuffle_full_solo(:,shuffle_i,bin_i)=diag(PV_corr_shuffle_full_solo);
                    %             corr_shuffle_full_solo(shuffle_i,bin_i)=nanmean(PV_corr_shuffle_full_solo(:,shuffle_i,bin_i));
                end
                
                
            end
            
            
            %% plots
            % figure prop:
            figure('units','normalized','outerposition',[0 0 1 1])
            %pos:
            width=0.08;
            height=0.08;
            hor_dis=0.1;
            ver_dis=0.13;
            
            x_pos(1)=0.1;
            x_pos(2)=x_pos(1)+hor_dis;
            x_pos(3)=x_pos(2)+hor_dis;
            x_pos(4)=x_pos(3)+hor_dis;
            x_pos(5)=x_pos(4)+hor_dis;
            x_pos(6)=x_pos(5)+hor_dis;
            x_pos(7)=x_pos(6)+hor_dis;
            x_pos(8)=x_pos(7)+hor_dis;
            x_pos(9)=x_pos(8)+hor_dis;
            
            y_pos(1)=0.8;
            y_pos(2)=y_pos(1)-ver_dis;
            y_pos(3)=y_pos(2)-ver_dis;
            y_pos(4)=y_pos(3)-ver_dis;
            y_pos(5)=y_pos(4)-ver_dis;
            y_pos(6)=y_pos(5)-ver_dis;
            
            %1. plot matrix of cells by position for solo
            axes('position',[x_pos(1),y_pos(1),width,height])
            normalize_to_other_max=[];
            X_bins_vector=1:size(PV_solo_dir,2);
            X_bins_vector_of_centers=0.5:0.5:size(PV_solo_dir,2)-0.5;
            [xlimits, ylimits] = fn_plot_2D_field (PV_solo_dir, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
            xlabel('cells')
            ylabel('Position in tunnel (m)')
            title('tuning curves of all cells')
            
            axes('position',[x_pos(2),y_pos(1),width,height])
            normalize_to_other_max=[];
            X_bins_vector=1:size(normMAX_PV_solo_dir,2);
            X_bins_vector_of_centers=0.5:0.5:size(normMAX_PV_solo_dir,2)-0.5;
            [xlimits, ylimits] = fn_plot_2D_field (normMAX_PV_solo_dir, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
            % xlabel('cells')
            % ylabel('Position in tunnel (m)')
            title('Normalized to max')
            
            axes('position',[x_pos(3),y_pos(1),width,height])
            normalize_to_other_max=[];
            X_bins_vector=1:size(normMEAN_PV_solo_dir,2);
            X_bins_vector_of_centers=0.5:0.5:size(normMEAN_PV_solo_dir,2)-0.5;
            [xlimits, ylimits] = fn_plot_2D_field (normMEAN_PV_solo_dir, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
            %  xlabel('cells')
            %  ylabel('Position in tunnel (m)')
            title('Normalized to mean')
            
            % 2. plot matrices of cells by position for co for all bins:
            relevant_bins=1:4:length(bin_start);
            for bin_i=1:length(relevant_bins)
                PV_co_dir_bin_i=PV_co_dir_bin(:,:,relevant_bins(bin_i));
                normMAX_PV_co_dir_bin_i=PV_co_dir_bin_i./(max(PV_co_dir_bin_i));
                axes('position',[x_pos(bin_i),y_pos(2),width,height])
                normalize_to_other_max=[];
                X_bins_vector=1:size(normMAX_PV_co_dir_bin_i,2);
                X_bins_vector_of_centers=0.5:0.5:size(normMAX_PV_co_dir_bin_i,2)-0.5;
                [xlimits, ylimits] = fn_plot_2D_field (normMAX_PV_co_dir_bin_i, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
                
                if bin_i~=1
                    xticks([])
                    yticks([])
                else
                    ylabel('co bin mat (norm to max)')
                end
                
                axes('position',[x_pos(bin_i),y_pos(3),width,height])
                shuffle_i=randi(n_shuffles);
                x=remove_data_tuning_curve(:,:,shuffle_i,relevant_bins(bin_i));
                norm_x=x./(max(x+eps));
                normalize_to_other_max=[];
                X_bins_vector=1:size(norm_x,2);
                X_bins_vector_of_centers=0.5:0.5:size(norm_x,2)-0.5;
                [xlimits, ylimits] = fn_plot_2D_field (norm_x, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
                if bin_i~=1
                    xticks([])
                    yticks([])
                else
                    ylabel('shuffle (solo) bin mat (norm to max)')
                    
                end
                %             axes('position',[x_pos(bin_i),y_pos(4),width,height])
                %             shuffle_i=randi(n_shuffles);
                %             x=rest_of_data_tuning_curve(:,:,shuffle_i,relevant_bins(bin_i));
                %             norm_x=x./(max(x+eps));
                %             normalize_to_other_max=[];
                %             X_bins_vector=1:size(norm_x,2);
                %             X_bins_vector_of_centers=0.5:0.5:size(norm_x,2)-0.5;
                %             [xlimits, ylimits] = fn_plot_2D_field (norm_x, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
                %             if bin_i~=1
                %                 xticks([])
                %                 yticks([])
                %                 else
                %                 ylabel('rest of solo bin mat (norm to max)')
                
                %end
            end
            
            % 3. plot matrices of PV corr for all bins:
            %         axes('position',[x_pos(1),y_pos(3),width,height])
            %         normalize_to_other_max=[];
            %         X_bins_vector=1:size(PV_corr_bin,2);
            %         X_bins_vector_of_centers=0.5:0.5:size(PV_corr_bin,2)-0.5;
            %         [xlimits, ylimits] = fn_plot_2D_field (PV_corr_bin, X_bins_vector, X_bins_vector_of_centers,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
            %         % 3. plot matrices of cell corr for all bins:
            %         axes('position',[x_pos(1),y_pos(6),width,height])
            %         normalize_to_other_max=[];
            %         X_bins_vector=1:size(cell_corr_bin,2);
            %         X_bins_vector_of_centers=0.5:0.5:size(cell_corr_bin,2)-0.5;
            %         cells_bin=1:size(cell_corr_bin,1);
            %         cells_bin_vector_of_centers=0.5:0.5:size(cell_corr_bin,1)-0.5;
            %         [xlimits, ylimits] = fn_plot_2D_field (cell_corr_bin, X_bins_vector, X_bins_vector_of_centers,cells_bin, cells_bin_vector_of_centers,normalize_to_other_max);
            %
            %  save:
            fig_name=fullfile(pop_vec_folder,['pop_vector_mats_dir_',num2str(dir_i),'_binsize_',num2str(solo_X_bin_size),'min_spike',num2str(min_n_spike),'intersect_data_',num2str(run_only_on_intersent_data),'full_intersect',num2str(full_data),'.png']);
            
            saveas(gcf,fig_name)
            %% second fig
            % figure prop:
            figure('units','normalized','outerposition',[0 0 1 1])
            %pos:
            width=.1;
            height=.1;
            hor_dis=0.2;
            ver_dis=0.18;
            
            x_pos(1)=0.08;
            x_pos(2)=x_pos(1)+hor_dis;
            x_pos(3)=x_pos(2)+hor_dis;
            x_pos(4)=x_pos(3)+hor_dis;
            x_pos(5)=x_pos(4)+hor_dis;
            x_pos(6)=x_pos(5)+hor_dis;
            x_pos(7)=x_pos(6)+hor_dis;
            x_pos(8)=x_pos(7)+hor_dis;
            
            y_pos(1)=0.8;
            y_pos(2)=y_pos(1)-ver_dis*0.8;
            y_pos(3)=y_pos(2)-ver_dis;
            y_pos(4)=y_pos(3)-ver_dis*0.8;
            y_pos(5)=y_pos(4)-ver_dis;
            y_pos(6)=y_pos(5)-ver_dis;
            
            %1. plot hist of cells corr
            %     x_vec=-1:0.1:1;
            %     axes('position',[x_pos(dir_i),y_pos(1),width,height])
            %     h=hist(nanmean(cell_corr_shuffle_i(:,4:5),2),x_vec);
            %     b1 = bar(x_vec,h/sum(h),'FaceColor',[0.5,0.5,0.5]);
            %     set(get(b1,'Children'),'FaceAlpha',0.3)
            %     hold on
            %     h=hist(nanmean(cell_corr_bin(:,4:5),2),x_vec);
            %     b2 = bar(x_vec,h/sum(h),'FaceColor',[255,0,255]./255);
            %     set(get(b2,'Children'),'FaceAlpha',0.2)
            %     title('Cell correlations')
            %     l=legend('shuffle (solo vs solo)','solo vs co (only in +-10)','location','eastoutside');
            %     set(l, 'position',[x_pos(dir_i)+1.05*width y_pos(1) 0.1 0.05])
            %     xlabel('r')
            %     ylabel('proportion of cells')
            %     % plot cell cor by bin vs shuffle
            %     axes('position',[x_pos(1),y_pos(2),width,height])
            %     plot(bin_center,cell_corr_mean_shuffle,'color',[.5 .5 .5])
            %     hold on;
            %     plot(bin_center,cell_corr_mean_bin,'linewidth',2,'color','m')
            %     ylabel('r')
            %     xlabel('Inter-bat distance (m)')
            %     title('cell corr')
            %     % plot cell corr by bin vs shuffle
            %     axes('position',[x_pos(1),y_pos(3),width,height])
            %     plot(bin_center,corr_shuffle,'color',[.5 .5 .5])
            %     hold on;
            %     plot(bin_center,corr_bin,'linewidth',2,'color','m')
            %     ylabel('r')
            %     xlabel('Inter-bat distance (m)')
            %     title('Population vector corr')
            ylim_corr=[-0.2 0.7];
            % shuffle vs full solo
            %noraml:
            % plot cell cor by bin vs shuffle
            axes('position',[x_pos(1),y_pos(1),width,height])
            plot(bin_center,cell_corr_mean_bin_shuffle,'color',[.5 .5 .5])
            hold on;
            plot(bin_center,cell_corr_mean_bin_shuffle_lower,'color',[.5 .5 .5])
            plot(bin_center,other_dir_cell_corr_mean_bin,'color',[0 0 0])
            plot(bin_center,cell_corr_mean_bin,'linewidth',2,'color','m')
            ylim(ylim_corr)
            ylabel('r')
            xlabel('Inter-bat distance (m)')
            title('corr of PSTH per cell averaged across cells')
            
            % plot cell corr mat
            axes('position',[x_pos(1),y_pos(2),width,height])
            normalize_to_other_max=[];
            cells_bin=1:size(cell_corr_bin,1);
            cells_bin_vector_of_centers=0.5:0.5:size(cell_corr_bin,1)-0.5;
            [xlimits, ylimits] = fn_plot_2D_field (cell_corr_bin, bin_dis, bin_center,cells_bin, cells_bin_vector_of_centers,normalize_to_other_max);
            xlabel('Inter-bat distance (m)')
            ylabel('cells')
            title('corr of PSTH per cell')
            
            % plot cell corr by bin vs shuffle
            axes('position',[x_pos(1),y_pos(3),width,height])
            plot(bin_center,PV_corr_mean_bin_shuffle,'color',[.5 .5 .5])
            hold on;
            plot(bin_center,PV_corr_mean_bin_shuffle_lower,'color',[.5 .5 .5])
            plot(bin_center,other_dir_PV_corr_mean_bin,'color',[0 0 0])
            plot(bin_center,PV_corr_mean_bin,'linewidth',2,'color','m')
            ylim(ylim_corr)
            ylabel('r')
            xlabel('Inter-bat distance (m)')
            title(sprintf('Corr of allo population vector\nper position averaged across positions\n- no normalization'))
            
            % plot PV corr mat
            axes('position',[x_pos(1),y_pos(4),width,height])
            normalize_to_other_max=[];
            X_bins_vector=1:size(PV_corr_bin,2);
            X_bins_vector_of_centers=0.5:0.5:size(PV_corr_bin,2)-0.5;
            [xlimits, ylimits] = fn_plot_2D_field (PV_corr_bin, bin_dis, bin_center, solo_X_bins_vector,solo_X_bins_vector_of_centers,normalize_to_other_max);
            xlabel('Inter-bat distance (m)')
            ylabel('Position in tunnel (m)')
            title('corr of population vector per position - no normaliztion')
            
            %norm to max:
            % plot cell corr by bin vs shuffle
            axes('position',[x_pos(2),y_pos(3),width,height])
            plot(bin_center,normMAX_PV_corr_mean_bin_shuffle,'color',[.5 .5 .5])
            hold on;
            plot(bin_center,normMAX_PV_corr_mean_bin_shuffle_lower,'color',[.5 .5 .5])
            plot(bin_center,other_dir_normMAX_PV_corr_mean_bin,'color',[0 0 0])
            plot(bin_center,normMAX_PV_corr_mean_bin,'linewidth',2,'color','m')
            ylim(ylim_corr)
            ylabel('r')
            xlabel('Inter-bat distance (m)')
            title('Population vector corr - norm to max')
            % plot PV corr mat
            axes('position',[x_pos(2),y_pos(4),width,height])
            normalize_to_other_max=[];
            X_bins_vector=1:size(normMAX_PV_corr_bin,2);
            X_bins_vector_of_centers=0.5:0.5:size(normMAX_PV_corr_bin,2)-0.5;
            [xlimits, ylimits] = fn_plot_2D_field (normMAX_PV_corr_bin, bin_dis, bin_center, solo_X_bins_vector,solo_X_bins_vector_of_centers,normalize_to_other_max);
            xlabel('Inter-bat distance (m)')
            ylabel('Position in tunnel (m)')
            title('corr of population vector per position - norm to max')
            
            %norm to mean:
            
            % plot cell corr by bin vs shuffle
            axes('position',[x_pos(3),y_pos(3),width,height])
            plot(bin_center,normMEAN_PV_corr_mean_bin_shuffle,'color',[.5 .5 .5])
            hold on;
            plot(bin_center,normMEAN_PV_corr_mean_bin_shuffle_lower,'color',[.5 .5 .5])
            plot(bin_center,other_dir_normMEAN_PV_corr_mean_bin,'color',[0 0 0])
            plot(bin_center,normMEAN_PV_corr_mean_bin,'linewidth',2,'color','m')
            ylim(ylim_corr)
            ylabel('r')
            xlabel('Inter-bat distance (m)')
            title('Population vector corr - norm to mean')
            
            % plot PV corr mat
            axes('position',[x_pos(3),y_pos(4),width,height])
            normalize_to_other_max=[];
            [xlimits, ylimits] = fn_plot_2D_field (normMEAN_PV_corr_bin, bin_dis, bin_center,solo_X_bins_vector, solo_X_bins_vector_of_centers,normalize_to_other_max);
            xlabel('Inter-bat distance (m)')
            ylabel('Position in tunnel (m)')
            title('corr of population vector per position - norm to mean')
            
            % plot SI
            axes('position',[x_pos(1),y_pos(5),width,height])
            mean_shuffle=squeeze(nanmean(information_per_spike_shuffle));
            plot(bin_center,mean_shuffle,'color',[.5 .5 .5],'linewidth',2)
            hold on;
            plot(bin_center,nanmean(information_per_spike_co),'linewidth',2,'color','m')
            ylim([0 max([max(mean_shuffle)';nanmean(information_per_spike_co)'])])
            ylabel('SI')
            xlabel('Inter-bat distance (m)')
            
            % plot velocity
            axes('position',[x_pos(2),y_pos(5),width,height])
            plot(bin_dis,nanmean(vel_by_dis),'linewidth',2,'color','m')
            ylim([0 max(nanmean(vel_by_dis))*1.05])
            xlim([min(bin_dis) max(bin_dis)])
            xlabel('Inter-bat distance (m)')
            ylabel('Velocity (m/s)')
            
            % plot y diff
            axes('position',[x_pos(3),y_pos(5),width,height])
            plot(bin_dis,nanmean(y_pos_by_dis),'linewidth',2,'color','m')
            ylim([min(nanmean(y_pos_by_dis))*1.05 max(nanmean(y_pos_by_dis))*1.05])
            xlim([min(bin_dis) max(bin_dis)])
            xlabel('Inter-bat distance (m)')
            ylabel('relative deviation in y')
            title('deviation from mean y pos in solo per x pos')
            
            % plot intercest/union 10-90%
            
            intersect_by_union(intersect_by_union==0)=nan;
            axes('position',[x_pos(4),y_pos(5),width,height])
            [h x]=hist(intersect_by_union)
            bar(x,h/sum(h))
            title(['10-90 intersect/union',num2str(nanmean(intersect_by_union(:)))])
            xlim([0 1])
            
            %full data
            full_intersect_by_union(full_intersect_by_union==0)=nan;
            axes('position',[x_pos(5),y_pos(5),width,height])
            [h x]=hist(full_intersect_by_union)
            bar(x,h/sum(h))
            title(['full data intersect/union',num2str(nanmean(full_intersect_by_union(:)))])
            xlim([0 1])
            
            %save figure
            fig_name=fullfile(pop_vec_folder,['correctionpop_vector_fig_dir_',num2str(dir_i),'_binsize_',num2str(solo_X_bin_size),'min_spike',num2str(min_n_spike),'intersect_data_',num2str(run_only_on_intersent_data),'full_intersect',num2str(full_data),'.png']);
            saveas(gcf,fig_name)
            
        end
    end
end
