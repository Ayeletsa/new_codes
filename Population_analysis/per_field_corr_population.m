% population per field
clear
% parameters:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis_k_1.5_th_1\';
dir_info=dir(dir_data);

co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct\';
figure_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\population\';
%
min_n_spike=100;
SI_threshold=1;
min_dis_pos_neg=0.25*1e6;%us
min_n_spike_per_field=50;
min_r_length_per_field=0.5; %min length of non nan firing rate in per field
per_field_tuning_curve_all_cells=[];
cell_count=0;
sig_bin_count=0;
fr_count_all=0;

%% params for valid cells during CO
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
param_file_name=fullfile(param_folder,'co_population_params.mat');
load(param_file_name)
ego_signif_cells=[];
per_field_signif=[];
n_smooth=3;
all_pos_width=cell(n_smooth,1);
all_neg_width=cell(n_smooth,1);
all_compound_width=cell(n_smooth,1);
all_compound_pos_width=cell(n_smooth,1);
all_compound_neg_width=cell(n_smooth,1);
all_pos_rise_time=cell(n_smooth,1);
all_neg_rise_time=cell(n_smooth,1);
%% TO DO:
%check why there are large fields
%hist of per field sizes
% set minimum spikes for per field analysis
%%
for cell_i=3:length({dir_info.name})-1
    for dir_i=1:2
        % load data:
        load(fullfile(dir_data,dir_info(cell_i).name))
        bat=cell_co_solo_initial_analysis.exp_data.bat;
        day=cell_co_solo_initial_analysis.exp_data.day;
        cell_num=cell_co_solo_initial_analysis.exp_data.cell_num;
        
        all_ego_tuning(cell_i-2,dir_i,:)=cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:);
        num_spike_during_flight=sum(~isnan(cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(dir_i).spikes.ts_usec(:)));
        %if num_spike_during_flight<=min_n_spike || cell_co_solo_initial_analysis.solo(dir_i).SI<SI_threshold
        %   continue
        %%end
        
        
        %% ego signif cell
        shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
        file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);
        
        load(file_name)
        
        n_spikes = cell_co_solo_initial_analysis.co(dir_i).info.n_spikes;
        even_odd_coherence = shuffling_struct(dir_i).odd_even_coherence.corr;
        ego_inf = shuffling_struct(dir_i).shuffled_data.params.information_per_spike_ego.values(1);
        ego_inf_p = shuffling_struct(dir_i).shuffled_data.params.information_per_spike_ego.p_value;
        cv_p=shuffling_struct(dir_i).shuffled_data.params.cv.p_value;
        
        % conditions :
        a = n_spikes > min_spikes;
        b = even_odd_coherence > min_even_odd;
        c = ego_inf > min_ego_inf;
        d = ego_inf_p < alpha;
        %d=cv_p<alpha;
        e=cell_co_solo_initial_analysis.exp_data.mean_fr<max_for_pyramidal;
        
        if min([a,b,c,d,e]) == 1
            ego_signif_cells(cell_i-2,dir_i)=1;
        else
            ego_signif_cells(cell_i-2,dir_i)=0;
        end
        
        
        %% per field data
        if  ~isempty([cell_co_solo_initial_analysis.solo(dir_i).field_height]) && e==1
            cell_count=cell_count+1;
            Isolation_dis(cell_count)=cell_co_solo_initial_analysis.exp_data.Isolation_dis;
            L_Ratio(cell_count)=cell_co_solo_initial_analysis.exp_data.L_Ratio;
            fr_count=0;
            for fr_i=1:length([cell_co_solo_initial_analysis.solo(dir_i).field_height])
                if sum(isnan(cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_x_fr_per_field{1, fr_i}))<0.5*length(cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_x_fr_per_field{1, fr_i})
                    fr_count=fr_count+1;
                    fr_count_all=fr_count_all+1;
                    per_field_tuning_curve_all_cells{cell_count,fr_count}=cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_x_fr_per_field{1, fr_i}';
                    firing_rate=cell_co_solo_initial_analysis.co(dir_i).firing_rate;
                    solo=cell_co_solo_initial_analysis.solo(dir_i);
                    per_field_SI(fr_count_all)=firing_rate.information_per_spike_per_field{fr_i};
                    per_field_CV(fr_count_all)=firing_rate.cv{fr_i};
                    per_field_modulation_depth(fr_count_all)=firing_rate.modulation_depth{fr_i};
                    per_field_sparsity(fr_count_all)=firing_rate.sparsity{fr_i};
                    per_field_n_spikes(fr_count_all)=firing_rate.number_of_spikes_per_field{fr_i};
                    per_field_height(fr_count_all)=solo.field_height(fr_i);
                    per_field_nrm_height(fr_count_all)=solo.field_height_norm_by_mean(fr_i);
                    per_field_size(fr_count_all)=solo.field_size(fr_i);
                    
                    
                    
                    
                    per_field_signif(cell_count)=0;% it will change if it will be true
                    for fr_i=1:length(firing_rate.time_signif_field_new_smooth)
                        
                        if iscell(firing_rate.time_signif_field_new_smooth{fr_i}) & firing_rate.number_of_spikes_per_field{fr_i}>min_n_spike_per_field
                            
                            for smooth_i=1:n_smooth
                                r=firing_rate.time_fr_per_field_new_smooth_smooth_window{fr_i}(smooth_i,:);
                                non_nan_ind_r=find(~isnan(r));
                                [ind_length,~,~]=find_length_of_consecutive_ind(non_nan_ind_r,length(r));
                                
                                if firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.signif_based_on_extreme_bins==1 & ind_length>min_r_length_per_field*length(r)
                                    
                                    per_field_signif(cell_count)=1;
                                    neg_width=firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.width_neg_interp;
                                    neg_width(isnan(neg_width))=[];
                                    pos_width=firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.width_pos_interp;
                                    pos_width(isnan(pos_width))=[];
                                    pos_rise_time=firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.pos_rise_time_interp;
                                    pos_rise_time(isnan(pos_rise_time))=[];
                                    neg_rise_time=firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.neg_rise_time_interp;
                                    neg_rise_time(isnan(neg_rise_time))=[];
                                    
                                    comb_width=[];
                                    pos_comb=[];
                                    neg_comb=[];
                                    
                                    if ~isempty(pos_width) & ~isempty(neg_width) % if there are both positive and negative
                                        x_pos=firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.width_line_x_pos_interp;
                                        x_neg=firing_rate.time_signif_field_new_smooth{fr_i}{smooth_i}.width_line_x_neg_interp;
                                        x_pos(isnan(x_pos))=[];x_neg(isnan(x_neg))=[];
                                        % check if there are clsoe to each
                                        % other:
                                        if abs(min(x_pos(:))-max(x_neg(:)))<=min_dis_pos_neg | abs(max(x_pos(:))-min(x_neg(:)))<=min_dis_pos_neg
                                            %TO DO: to correct it to work
                                            %properly on all type of tuning
                                            %curve if we will use it!
                                            comb_width=[max([x_pos(:);x_neg(:)])-min([x_pos(:);x_neg(:)])];
                                            pos_comb=pos_width;
                                            neg_comb=neg_width;
                                            pos_width=[];
                                            neg_width=[];
                                            neg_rise_time=[];
                                            pos_rise_time=[];
                                        end
                                    end
                                    all_pos_width{smooth_i}=[all_pos_width{smooth_i},pos_width];
                                    all_neg_width{smooth_i}=[all_neg_width{smooth_i},neg_width];
                                    all_compound_width{smooth_i}=[all_compound_width{smooth_i},comb_width];
                                    all_compound_pos_width{smooth_i}=[all_compound_pos_width{smooth_i},pos_comb];
                                    all_compound_neg_width{smooth_i}=[all_compound_neg_width{smooth_i},neg_comb];
                                    all_pos_rise_time{smooth_i}=[all_pos_rise_time{smooth_i},pos_rise_time];
                                    all_neg_rise_time{smooth_i}=[all_neg_rise_time{smooth_i},neg_rise_time];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%%
% take only cells with more than one field to calculat:
[r c]=find(~cellfun(@isempty,per_field_tuning_curve_all_cells));
ind=zeros(size(per_field_tuning_curve_all_cells));
ind(sub2ind(size(per_field_tuning_curve_all_cells),r,c))=1;
cells_to_remove=sum(ind,2)<=1;
L_Ratio_relevant_cells=L_Ratio(~cells_to_remove);
Isolation_dis_relevant_cells=Isolation_dis(~cells_to_remove);
per_field_tuning_curve_relevant_cells=per_field_tuning_curve_all_cells(~cells_to_remove,:);
all_cells_corr=[];
other_cells_corr=[];
for cell_i=1:length(per_field_tuning_curve_relevant_cells)
    cell_vec=zeros(length(per_field_tuning_curve_relevant_cells),1);
    cell_vec(cell_i)=1;
    mat_1=[per_field_tuning_curve_relevant_cells{cell_i,:}];
    cell_corr_values=corr(mat_1,mat_1,'rows','pairwise');
    cell_corr_values=triu(cell_corr_values,1);
    cell_corr_values=cell_corr_values(find(cell_corr_values~=0));
    all_cells_corr=[all_cells_corr; cell_corr_values];
    mean_corr_cell_value(cell_i)=nanmean(cell_corr_values);
    
    other_cells_mat=[per_field_tuning_curve_relevant_cells{~cell_vec,:}];
    
    other_cell_corr_values=corr(mat_1,other_cells_mat,'rows','pairwise');
    other_cells_corr=[other_cells_corr;other_cell_corr_values(:)];
end

% only signif per field:
relevant_cells=~cells_to_remove & per_field_signif';
per_field_tuning_curve_relevant_cells=per_field_tuning_curve_all_cells(relevant_cells,:);
cells_with_signif_bin_corr=[];
other_cells_with_signif_bin_cells_corr=[];
for cell_i=1:length(per_field_tuning_curve_relevant_cells)
    cell_vec=zeros(length(per_field_tuning_curve_relevant_cells),1);
    cell_vec(cell_i)=1;
    mat_1=[per_field_tuning_curve_relevant_cells{cell_i,:}];
    cell_corr_values=corr(mat_1,mat_1,'rows','pairwise');
    cell_corr_values=triu(cell_corr_values,1);
    cells_with_signif_bin_corr=[cells_with_signif_bin_corr; cell_corr_values(find(cell_corr_values~=0))];
    
    other_cells_mat=[per_field_tuning_curve_relevant_cells{~cell_vec,:}];
    
    other_cell_corr_values=corr(mat_1,other_cells_mat,'rows','pairwise');
    other_cells_with_signif_bin_cells_corr=[other_cells_with_signif_bin_cells_corr;other_cell_corr_values(:)];
end

% egocentric signif cells tuning
signif_ego_corr_cell=[];
signif_ego_corr_other_cells=[];
any_dir_signif_cells=sum(ego_signif_cells,2);
ego_signif_tuning=all_ego_tuning(any_dir_signif_cells>0,:,:);
for cell_i=1:size(ego_signif_tuning,1)
    
    mat_1=squeeze(ego_signif_tuning(cell_i,:,:));
    cell_corr_values=corr(mat_1(1,:)',mat_1(2,:)','rows','pairwise');
    signif_ego_corr_cell=[signif_ego_corr_cell,cell_corr_values];
    
    cell_vec=zeros(size(ego_signif_tuning,1),1);
    cell_vec(cell_i)=1;
    other_cells=ego_signif_tuning(~cell_vec,:,:);
    other_cells=[squeeze(other_cells(:,1,:));squeeze(other_cells(:,2,:))];
    other_cell_corr_values=corr(mat_1(1,:)',other_cells','rows','pairwise');
    signif_ego_corr_other_cells=[signif_ego_corr_other_cells;other_cell_corr_values(:)];
    other_cell_corr_values=corr(mat_1(2,:)',other_cells','rows','pairwise');
    signif_ego_corr_other_cells=[signif_ego_corr_other_cells;other_cell_corr_values(:)];
end

% all egocentric tuning
all_ego_corr_cell=[];
all_ego_corr_other_cells=[];
for cell_i=1:length(all_ego_tuning)
    mat_1=squeeze(all_ego_tuning(cell_i,:,:));
    cell_corr_values=corr(mat_1(1,:)',mat_1(2,:)','rows','pairwise');
    all_ego_corr_cell=[all_ego_corr_cell,cell_corr_values];
    
    cell_vec=zeros(length(all_ego_tuning),1);
    cell_vec(cell_i)=1;
    other_cells=all_ego_tuning(~cell_vec,:,:);
    other_cells=[squeeze(other_cells(:,1,:));squeeze(other_cells(:,2,:))];
    other_cell_corr_values=corr(mat_1(1,:)',other_cells','rows','pairwise');
    all_ego_corr_other_cells=[all_ego_corr_other_cells;other_cell_corr_values(:)];
    other_cell_corr_values=corr(mat_1(2,:)',other_cells','rows','pairwise');
    all_ego_corr_other_cells=[all_ego_corr_other_cells;other_cell_corr_values(:)];
end
%% PLOT

figure('units','normalized','outerposition',[0 0 1 1])

subplot(3,2,1)
plot_hists(all_ego_corr_cell,all_ego_corr_other_cells,'all cells egocentric corr (between directions)')

subplot(3,2,2)
plot_hists(signif_ego_corr_cell,signif_ego_corr_other_cells,'ego signif cells egocentric corr (between directions)')

subplot(3,2,3)
plot_hists(all_cells_corr,other_cells_corr,'per field egocentric correlations')

subplot(3,2,4)
plot_hists(cells_with_signif_bin_corr,other_cells_with_signif_bin_cells_corr,'per field egocentric correlations- cells with signif bin in at least one field')

subplot(3,2,5)
plot_scatter(mean_corr_cell_value',L_Ratio_relevant_cells','mean corr per field','L Ratio')

subplot(3,2,6)
plot_scatter(mean_corr_cell_value',Isolation_dis_relevant_cells','mean corr per field','Isolation distance')

% save figure:
file_name=fullfile(figure_folder,'hist_corr_directions_and_per_field.jpg');
saveas(gcf,file_name)
%% Figure of per field tuning with
figure('units','normalized','outerposition',[0 0 1 1])
smooth_vec=[3,5,7];

n_plots=7;
for smooth_i=1:n_smooth
    x_vec=(0:0.2:4);
    xlimits=[0 4];
    smooth_wind=smooth_vec(smooth_i);
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+1)
    data=all_pos_width{smooth_i}./1e6;
    txt='Positive fields:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+2)
    data=all_neg_width{smooth_i}./1e6;
    txt='Negative fields:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+3)
    data=all_compound_width{smooth_i}./1e6;
    txt='Compound fields:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+4)
    data=all_compound_pos_width{smooth_i}./1e6;
    txt='Pos from conpound:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+5)
    data=all_compound_neg_width{smooth_i}./1e6;
    txt='Neg from compound:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
    x_vec=(0:0.1:2);
    xlimits=[0 2];
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+6)
    data=all_pos_rise_time{smooth_i}./1e6;
    txt='Pos rise time:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
    subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+7)
    data=all_neg_rise_time{smooth_i}./1e6;
    txt='Neg rise time:';
    hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
    
end

file_name=fullfile(figure_folder,'hist_per_field_tuning_width.jpg');
saveas(gcf,file_name)
%%

figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,4,1)
plot_scatter(per_field_n_spikes',per_field_SI','# spikes','SI')
subplot(4,4,2)
plot_scatter(per_field_n_spikes',per_field_CV','# spikes','CV')
subplot(4,4,3)
plot_scatter(per_field_n_spikes',per_field_modulation_depth','# spikes','modulation depth')
subplot(4,4,4)
plot_scatter(per_field_n_spikes',per_field_sparsity','# spikes','sparsity')

subplot(4,4,5)
plot_scatter(per_field_height',per_field_SI','field height','SI')
subplot(4,4,6)
plot_scatter(per_field_height',per_field_CV','field height','CV')
subplot(4,4,7)
plot_scatter(per_field_height',per_field_modulation_depth','field height','modulation depth')
subplot(4,4,8)
plot_scatter(per_field_height',per_field_sparsity','field height','sparsity')

subplot(4,4,9)
plot_scatter(per_field_nrm_height',per_field_SI','norm field height','SI')
subplot(4,4,10)
plot_scatter(per_field_nrm_height',per_field_CV','norm field height','CV')
subplot(4,4,11)
plot_scatter(per_field_nrm_height',per_field_modulation_depth','norm field height','modulation depth')
subplot(4,4,12)
plot_scatter(per_field_nrm_height',per_field_sparsity','norm field height','sparsity')


subplot(4,4,13)
plot_scatter(per_field_size',per_field_SI','field size','SI')
subplot(4,4,14)
plot_scatter(per_field_size',per_field_CV','field size','CV')
subplot(4,4,15)
plot_scatter(per_field_size',per_field_modulation_depth','field size','modulation depth')
subplot(4,4,16)
plot_scatter(per_field_size',per_field_sparsity','field size','sparsity')

% save figure:
file_name=fullfile(figure_folder,'per_field_properties_corr.jpg');
saveas(gcf,file_name)


%% functions
function plot_scatter(x,y,x_name,y_name)
[r,p]=corr(x,y,'rows','pairwise');
scatter(x,y)
xlabel(x_name)
ylabel(y_name)
title(sprintf('r=%.2f p=%.2f',r,p))
end

function plot_hists(x_data,shuffle,name)
x_bin=-1:0.1:1;
[h,x]=hist(shuffle,x_bin);
b1=bar(x,h/sum(h));
b1.FaceAlpha = 0.5;
hold on
[h,x]=hist(x_data,x_bin);
b2=bar(x,h/sum(h));
b2.FaceAlpha = 0.5;
xlabel('r')
ylabel('proportion')

legend('shuffle','within a cell')
[h,p] = kstest2(x_data,shuffle);

title(sprintf('%s, p=%.2f', name,p))
end

function hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
[h,x]=hist(data,x_vec);
bar(x,h/sum(h))
n_fields=length(data);
xlabel('tuning width (s)')
ylabel('proportion')
mean_hist=nanmean(data);
median_hist=nanmedian(data);
title(sprintf('%s \nsmooth=%d #fields=%d\nmean=%.2f meadian=%.2f',txt,smooth_wind,n_fields,mean_hist,median_hist))
xlim(xlimits)
end