%% Population per field
clear

%% data:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_info=dir(dir_data);
analysis_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\population\';


co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct\';
figure_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\population\';

inclusion_dir='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\inclusion_cells_struct\';
inclusion_dir_info=dir(inclusion_dir);
inclusion_names={inclusion_dir_info.name};
inclusion_names([inclusion_dir_info.isdir])=[];
cell_ind=regexp(inclusion_names{2},'cell_');

cell_nums_inclusion=(cellfun(@(c) c(cell_ind+5:cell_ind+7),inclusion_names,'UniformOutput',false));

%% params for valid cells during CO
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
pop_param_file_name=fullfile(param_folder,'co_population_params.mat');
load(pop_param_file_name)
per_field_param_file_name=fullfile(param_folder,'per_field_params.mat');
load(per_field_param_file_name)
co_param_file_name=fullfile(param_folder,'co_params.mat');
load(co_param_file_name)

%% initialize:
per_field_cell_dir_count=0;
fr_count_all=0;
cell_dir_per_field_count=0;
fr_count_signif=0;

all_per_field_valid_tuning={};
all_pos_width=[];
all_neg_width=[];
all_compound_width=[];
all_compound_pos_width=[];
all_compound_neg_width=[];
all_pos_rise_time=[];
all_neg_rise_time=[];
all_compound_pos_rise_time=[];
all_compound_neg_rise_time=[];
all_peak_fr_pos_sig_bins=[];
all_ego_position_pos_sig_bins=[];
all_ego_position_neg_sig_bins=[];
%% TO DO:
%check why there are large fields
%hist of per field sizes
% set minimum spikes for per field analysis
%%
file_names={dir_info.name};
file_names=file_names(find([dir_info.isdir]==0));

%% 1. run over cells:
for cell_i=1:length(file_names)
    % load data:
    load(fullfile(dir_data,file_names{cell_i}))
    cell_num=cell_co_solo_initial_analysis.exp_data.cell_num;
    %load cell's inclusion:
    Match=cellfun(@(a) find(contains(a,num2str(cell_num))) , cell_nums_inclusion, 'UniformOutput', 0);
    r=find(cellfun(@(c) ~isempty(c),Match));
    load(fullfile(inclusion_dir,inclusion_names{r}))
    
    
    for dir_i=1:2
        
       % basic cell data:
        bat=cell_co_solo_initial_analysis.exp_data.bat;
        day=cell_co_solo_initial_analysis.exp_data.day;
        cell_num=cell_co_solo_initial_analysis.exp_data.cell_num;
        %ego tuning of the cell:
        all_ego_tuning(cell_i,dir_i,:)=cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:);
        all_zscore_ego_tuning(cell_i,dir_i,:)=zscore(cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:));

        %% ego signif cell
 
%         shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
%         shuffle_file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);

        ego_signif_cells(cell_i,dir_i)=inclusion(dir_i).ego_cell;
        pyramidal_cells(cell_i,dir_i)=inclusion(dir_i).pyr;
        %% per field data
        %1. Test if run on per field:
       
%         a=sum(~isnan(cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(dir_i).spikes.ts_usec(:)))>=min_n_spike;
%         b=cell_co_solo_initial_analysis.solo(dir_i).SI>SI_threshold;
%         c=~isempty(cell_co_solo_initial_analysis.solo(dir_i).fields);
%         d=cell_co_solo_initial_analysis.exp_data.mean_fr<max_for_pyramidal;
%         per_field_cell_cond=a & b & c & d;
        if inclusion(dir_i).place_cell 
          cell_dir_per_field_count=cell_dir_per_field_count+1;

            per_field_data=cell_co_solo_initial_analysis.co(dir_i).per_field_href;
            solo_field_data=cell_co_solo_initial_analysis.solo(dir_i).fields;
            
          

          for field_i=1:length(cell_co_solo_initial_analysis.co(dir_i).per_field_href)
              % run only on per fields that obey conditions:
              r=per_field_data(field_i).tuning_dis_x_fr_per_field;
              [ind_length,~,~]=find_length_of_consecutive_ind(find(~isnan(r)),length(r));
              a=max(ind_length)>=min_r_length_per_field*length(r); %enough data that is not nan in tuning curve
              b=per_field_data(field_i).number_of_spikes_per_field>=min_n_spike_per_field;% enough spikes
              per_field_cond=a & b;
               if per_field_cond
                fr_count_all=fr_count_all+1;
                %1. per field tuning:
                all_per_field_valid_tuning{cell_dir_per_field_count,field_i}=r';
                % iso dis:
                Isolation_dis(cell_dir_per_field_count)=cell_co_solo_initial_analysis.exp_data.Isolation_dis;
                % n fields:
                n_fields_per_cell_for_cell_scatter(cell_dir_per_field_count)=length([cell_co_solo_initial_analysis.solo(dir_i).fields.loc]);

                n_fields_per_cell_for_field_scatter(fr_count_all)=length([cell_co_solo_initial_analysis.solo(dir_i).fields.loc]);
                
                all_per_field_peak_solo_pos{cell_dir_per_field_count,field_i}=solo_field_data(field_i).loc;
                all_per_field_peak_solo{cell_dir_per_field_count,field_i}=solo_field_data(field_i).peak;
                all_per_field_width_solo{cell_dir_per_field_count,field_i}=solo_field_data(field_i).width_href;

                %2. get solo prop per field:
                per_field_solo_height(fr_count_all)=solo_field_data(field_i).peak;
                per_field_solo_nrm_height(fr_count_all)=solo_field_data(field_i).peak./nanmean(cell_co_solo_initial_analysis.solo(dir_i).PSTH_for_field_detection);
                per_field_solo_size(fr_count_all)=solo_field_data(field_i).width_href;
                per_field_solo_SI(fr_count_all)=solo_field_data(field_i).field_SI;
              
                zscore_r=zscore(r);
                %3. get co prop per field:
                per_field_co_SI(fr_count_all)=per_field_data(field_i).SI;
                per_field_co_CV(fr_count_all)=nanstd(r)./nanmean(r);
                per_field_co_modulation_depth(fr_count_all)=(max(r)-nanmean(r))./max(r);
                per_field_co_modulation_depth_zscore_r(fr_count_all)=(max(zscore_r)-min(zscore_r))./max(zscore_r);
                %per_field_co_modulation_depth_zscore_r(fr_count_all)=(max(zscore_r)-min(zscore_r));
                per_field_co_neg_modulation_depth(fr_count_all)=(nanmean(r)-min(r))./nanmean(r);
                per_field_co_iqr_by_median(fr_count_all)=iqr(r)./nanmedian(r);
                per_field_co_sparsity(fr_count_all)=per_field_data(field_i).sparsity;
                
                
                %4. get tuning width and rise time:
                if per_field_data(field_i).time_signif_field.signif_based_on_extreme_bins==1 
                    rise_and_width_data=get_per_field_tuning_width_and_rise_time(per_field_data(field_i).time_signif_field,min_dis_pos_neg);
                    id_per_field_non_signif(fr_count_all)=0;
                    field_signif_per_cell{cell_dir_per_field_count,field_i}=1;
                    id_per_field_pos_tuning(fr_count_all)=rise_and_width_data.pos_per_field_tuning;
                    id_per_field_neg_tuning(fr_count_all)=rise_and_width_data.neg_per_field_tuning;
                    

                    all_pos_width=[all_pos_width,rise_and_width_data.pos_width];
                    all_neg_width=[all_neg_width,rise_and_width_data.neg_width];
                    all_compound_width=[all_compound_width,rise_and_width_data.compound_width];
                    all_compound_pos_width=[all_compound_pos_width,rise_and_width_data.pos_compound_width];
                    all_compound_neg_width=[all_compound_neg_width,rise_and_width_data.neg_compound_width];
                    all_pos_rise_time=[all_pos_rise_time,rise_and_width_data.pos_rise_time];
                    all_neg_rise_time=[all_neg_rise_time,rise_and_width_data.neg_rise_time];
                    all_compound_pos_rise_time=[all_compound_pos_rise_time,rise_and_width_data.pos_rise_time_compund];
                    all_compound_neg_rise_time=[all_compound_neg_rise_time,rise_and_width_data.neg_rise_time_compund];
                    if per_field_data(field_i).dis_signif_field.signif_based_on_extreme_bins==1  
                    all_peak_fr_pos_sig_bins=[all_peak_fr_pos_sig_bins,per_field_data(field_i).tuning_dis_x_fr_per_field([per_field_data(field_i).dis_signif_field.pos_signif])];
                    all_ego_position_pos_sig_bins=[all_ego_position_pos_sig_bins,dis_per_field_bin_vec_of_center([per_field_data(field_i).dis_signif_field.pos_signif])]; 
                    all_ego_position_neg_sig_bins=[all_ego_position_neg_sig_bins,dis_per_field_bin_vec_of_center([per_field_data(field_i).dis_signif_field.neg_signif])]; 
                    %
                    fr_count_signif=fr_count_signif+1;
                    pos_signif_position_per_field{fr_count_signif}=dis_per_field_bin_vec_of_center([per_field_data(field_i).dis_signif_field.pos_signif]);
                    neg_signif_position_per_field{fr_count_signif}=dis_per_field_bin_vec_of_center([per_field_data(field_i).dis_signif_field.neg_signif]);
                    pos_signif_zscore_per_field{fr_count_signif}=zscore_r([per_field_data(field_i).dis_signif_field.pos_signif]);
                    neg_signif_zscore_per_field{fr_count_signif}=zscore_r([per_field_data(field_i).dis_signif_field.neg_signif]);
                    solo_field_hight_per_field_for_signif_fields(fr_count_signif)=solo_field_data(field_i).peak;
                    end
                    
                else
                    id_per_field_non_signif(fr_count_all)=1;
                    id_per_field_pos_tuning(fr_count_all)=0;
                    id_per_field_neg_tuning(fr_count_all)=0;
                    field_signif_per_cell{cell_dir_per_field_count,field_i}=0;
                end
               end
          end
%           % compute cell per field corr:
%           if size(all_per_field_valid_tuning,1)>1
%               per_field_cell_dir_count=per_field_cell_dir_count+1;
%               corrs=tril(corr(all_per_field_valid_tuning',all_per_field_valid_tuning','rows','pairwise'),-1);
%               corr_per_field_per_cell(per_field_cell_dir_count)=nanmean(corrs(corrs~=0));
%           
%           end
        end
    end
%     %% compute corr between dirs:
%     corrs=tril(corr(squeeze(all_ego_tuning(cell_i,:,:))','rows','pairwise'),-1);
%     corr_ego_tuning_between_dirs(cell_i)=corrs(2,1);
end
pos_signif_position_per_field_mat=convert_to_mat(pos_signif_position_per_field);
neg_signif_position_per_field_mat=convert_to_mat(neg_signif_position_per_field);
pos_signif_zscore_per_field_mat=convert_to_mat(pos_signif_zscore_per_field);
neg_signif_zscore_per_field_mat=convert_to_mat(neg_signif_zscore_per_field);
solo_field_hight_per_field_for_signif_fields_repmat_pos=repmat(solo_field_hight_per_field_for_signif_fields',1,size(pos_signif_position_per_field_mat,2));
solo_field_hight_per_field_for_signif_fields_repmat_neg=repmat(solo_field_hight_per_field_for_signif_fields',1,size(neg_signif_position_per_field_mat,2));

% figure; 
% subplot(1,2,1);scatter(solo_field_hight_per_field_for_signif_fields_repmat_pos(:),pos_signif_position_per_field_mat(:),[],pos_signif_zscore_per_field_mat(:));colorbar
% subplot(1,2,2);scatter(solo_field_hight_per_field_for_signif_fields_repmat_neg(:),neg_signif_position_per_field_mat(:),[],neg_signif_zscore_per_field_mat(:));colorbar

%% compute corrs and shuffles:

%1. cell ego tuning between directions:
%-------------------------------------- 
data=all_ego_tuning;
[all_cells_ego_corr_between,all_cells_ego_corr_within]=compute_corrs_within_and_between_cells(data);

%2. cell ego tuning between directions - only ego signif cells:
%--------------------------------------------------------------
%find signif cells:
ego_signif_cells_both_dir=sum(ego_signif_cells,2)>0;
data=all_ego_tuning(ego_signif_cells_both_dir,:,:);

[ego_signif_cells_ego_corr_between,ego_signif_cells_ego_corr_within]=compute_corrs_within_and_between_cells(data);
%3. cell per field ego tuning between directions - valid cells and fields:
%-------------------------------------------------------------------------
% a. find valid cells:
% take only cells with more than one field to calculate:
[r c]=find(~cellfun(@isempty,all_per_field_valid_tuning));
ind=zeros(size(all_per_field_valid_tuning));
ind(sub2ind(size(all_per_field_valid_tuning),r,c))=1;
cells_to_remove=sum(ind,2)<=1;
Isolation_dis_relevant_cells=Isolation_dis(~cells_to_remove);
n_fields_per_cell_relevant_cells=n_fields_per_cell_for_cell_scatter(~cells_to_remove);

per_field_tuning_curve_relevant_cells=all_per_field_valid_tuning(~cells_to_remove,:);
all_per_field_peak_solo_pos_relevant_cells=all_per_field_peak_solo_pos(~cells_to_remove,:);
all_per_field_peak_solo_relevant_cells=all_per_field_peak_solo(~cells_to_remove,:);
all_per_field_width_solo_relevant_cells =all_per_field_width_solo(~cells_to_remove,:);
% b. compute correlations within a cell and between cells:
data=per_field_tuning_curve_relevant_cells;
field_signif_per_cell=field_signif_per_cell(~cells_to_remove,:);
[per_field_cells_corr_between,per_field_cells_corr_within,per_field_mean_corr_within,dist_peak_solo,diff_peak_solo]=compute_corrs_within_and_between_cells_per_field(data,field_signif_per_cell,all_per_field_peak_solo_pos_relevant_cells,all_per_field_peak_solo_relevant_cells);
% [r,c]=find(cellfun(@(x) x==0,field_signif_per_cell));
% data([r,c])={};
% data=reshape(data,size(per_field_tuning_curve_relevant_cells,1),[])
n_fields_per_cell=sum(cellfun(@(a) ~isempty(a),data),2);
n_pairs=factorial(n_fields_per_cell)./(factorial(2).*factorial(n_fields_per_cell-2));
n_fields_per_cell_relevant_fields=repelem(n_fields_per_cell_relevant_cells,n_pairs);
isolation_dis_cell_relevant_fields=repelem(Isolation_dis_relevant_cells,n_pairs);

%% PLOT
%1. plot direction tuning corrs:
figure('units','normalized','outerposition',[0 0 1 1])
% ego tuning:
%-----------------
% corr between directions all cells
% subplot(4,7,1)
% plot_hists(all_cells_ego_corr_within,all_cells_ego_corr_between,'all cells egocentric corr (between directions)')
% 
% % corr between directions ego signif cells
% subplot(4,7,2)
% plot_hists(ego_signif_cells_ego_corr_within,ego_signif_cells_ego_corr_between,'ego signif cells egocentric corr (between directions)')
subplot(4,7,3)
for bin_i=1:length(dis_per_field_bin_vec_of_center)
    idx=find(all_ego_position_pos_sig_bins==dis_per_field_bin_vec_of_center(bin_i));
    mean_fr_per_bin(bin_i)=nanmean(all_peak_fr_pos_sig_bins(idx));
end
plot(dis_per_field_bin_vec_of_center,mean_fr_per_bin)
xlim([-40 40])
xlabel('Inter-bat distance (m)')
ylabel('mean FR')
title('mean FR for positive significant bins')

subplot(4,7,4)
x=all_ego_position_pos_sig_bins;
y=all_peak_fr_pos_sig_bins;
x_name='Inter-bat distance (m)';
y_name='Firing rate';
plot_scatter(x,y,x_name,y_name)
title(sprintf('position and firing rate\n of position significant bins'))
xlim([-40 40])

subplot(4,7,5)
h=hist(all_ego_position_pos_sig_bins,dis_per_field_bin_vec_of_center,'g');
bar(dis_per_field_bin_vec_of_center,h/sum(h))
xlabel('Inter-bat distance (m)')
ylabel('proportion')
title('pos signif bins')

subplot(4,7,6)
h=hist(all_ego_position_neg_sig_bins,dis_per_field_bin_vec_of_center,'r');
bar(dis_per_field_bin_vec_of_center,h/sum(h))
xlabel('Inter-bat distance (m)')
ylabel('proportion')
title('neg signif bins')

% [c,n]=histc(all_ego_position_pos_sig_bins,dis_per_field_bin_vec_of_center);
% results = cell(length(dis_per_field_bin_vec_of_center),1);
% for K = 1 : length(dis_per_field_bin_vec_of_center)
%   results{K} = all_peak_fr_pos_sig_bins( n == K );
% end
% mean_fr=cellfun(@mean,results);
% subplot(4,7,5)
% plot(dis_per_field_bin_vec_of_center,mean_fr)

subplot(4,7,7)
x=dist_peak_solo;
y=per_field_cells_corr_within;
x_name='distance peak solo';
y_name='per field corr';
plot_scatter(x,y,x_name,y_name)

subplot(4,7,14)
x=dist_peak_solo;
y=diff_peak_solo;
x_name='distance peak solo';
y_name='diff peak solo';
plot_scatter(x,y,x_name,y_name)

% subplot(4,7,7)
% x=diff_peak_solo;
% y=per_field_cells_corr_within;
% x_name='diff peak solo';
% y_name='per field corr';
% plot_scatter(x,y,x_name,y_name)

% % tuning ego cells
% subplot(4,7,3)
% pyramidal_cells_and_ego=ego_signif_cells&pyramidal_cells;
% 
% mat=[squeeze(all_zscore_ego_tuning(pyramidal_cells_and_ego(:,1),1,:));squeeze(all_zscore_ego_tuning(pyramidal_cells_and_ego(:,2),2,:))];
% [~,I]=max(mat,[],2);
% [~,I]=sort(I,'descend');
% mat=mat(I,:);
% imagesc(dis_X_bins_vector_of_centers,1:length(mat),mat)
% title('zscore ego signif & pyramidal cells')
% ylabel('cells')
% xlabel('Inter-bat distance (m)')
% 
% subplot(4,7,4)
% mat=[squeeze(all_zscore_ego_tuning(ego_signif_cells(:,1),1,:));squeeze(all_zscore_ego_tuning(ego_signif_cells(:,2),2,:))];
% [~,I]=max(mat,[],2);
% [~,I]=sort(I,'descend');
% mat=mat(I,:);
% imagesc(dis_X_bins_vector_of_centers,1:length(mat),mat)
% title('zscore ego signif cells')
% ylabel('cells')
% xlabel('Inter-bat distance (m)')
% colormap('jet')
% % plot mean
% % subplot(4,7,4)
% % plot(mean(mat))
% % tuning pyramidal ego cells
% subplot(4,7,5)
% non_pyramidal_cells_and_ego=ego_signif_cells&~pyramidal_cells;
% %non_pyramidal_cells_and_ego=~pyramidal_cells;
% mat=[squeeze(all_zscore_ego_tuning(non_pyramidal_cells_and_ego(:,1),1,:));squeeze(all_zscore_ego_tuning(non_pyramidal_cells_and_ego(:,2),2,:))];
% [~,I]=max(mat,[],2);
% [~,I]=sort(I,'descend');
% mat=mat(I,:);
% mat(nansum(mat,2)==0,:)=[];
% %imagesc(mat)
% imagesc(dis_X_bins_vector_of_centers,1:length(mat),mat)
% title('ego interneurons')
% 
% subplot(4,7,6)
% non_pyramidal=~pyramidal_cells;
% %non_pyramidal_cells_and_ego=~pyramidal_cells;
% mat=[squeeze(all_zscore_ego_tuning(non_pyramidal(:,1),1,:));squeeze(all_zscore_ego_tuning(non_pyramidal(:,2),2,:))];
% [~,I]=max(mat,[],2);
% [~,I]=sort(I,'descend');
% mat=mat(I,:);
% mat(nansum(mat,2)==0,:)=[];
% %imagesc(mat)
% imagesc(dis_X_bins_vector_of_centers,1:length(mat),mat)
% title('interneurons')

% Per field
%--------------------------------
% corr between per field for valid cells
subplot(4,7,8:9)
plot_hists(per_field_cells_corr_within,per_field_cells_corr_between,'per field egocentric correlations')

% scatter of isolation distance and per field mean corr:    
%subplot(4,7,10:11)
%plot_scatter(isolation_dis_cell_relevant_fields',per_field_cells_corr_within,'Isolation distance','pairs corr per field')

%scatter of n fields and mean corr
%subplot(4,7,12:13)
%plot_scatter(n_fields_per_cell_relevant_fields',per_field_cells_corr_within,'n fields','pairs corr per field')

% Tuning width and tise time -  per field:
%-----------------------------------------
x_vec=(0:0.2:4);
xlimits=[0 4];

subplot(4,7,15)
data=all_pos_width./1e6;
txt='Positive fields:';
x_str='Tuning width (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)
    
subplot(4,7,16)
data=all_neg_width./1e6;
txt='Negative fields:';
x_str='Tuning width (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)
    
subplot(4,7,17)
data=all_compound_width./1e6;
txt='Compound fields:';
x_str='Tuning width (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)
    
    
subplot(4,7,18)
data=all_compound_pos_width./1e6;
txt='Pos from conpound:';
x_str='Tuning width (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)
    
subplot(4,7,19)
data=all_compound_neg_width./1e6;
txt='Neg from compound:';
x_str='Tuning width (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)
    
x_vec=(0:0.1:2);
xlimits=[0 2];
    
subplot(4,7,20)
data=all_pos_rise_time./1e6;
txt='Pos rise time:';
x_str='rise time (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)
    
subplot(4,7,21)
data=all_neg_rise_time./1e6;
txt='Neg rise time:';
x_str='rise time (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)

subplot(4,7,27)
data=all_compound_pos_rise_time./1e6;
txt='pos from combine rise time:';
x_str='rise time (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)

subplot(4,7,28)
data=all_compound_neg_rise_time./1e6;
txt='neg from combine rise time:';
x_str='rise time (s)';
hist_per_field(data,x_vec,txt,xlimits,x_str)


% save figure:
file_name=fullfile(figure_folder,'hist_corr_directions_and_per_field.jpg');
saveas(gcf,file_name)

%% plot scatters of per field


per_field_co_SI=real(per_field_co_SI);
per_field_co_SI(find(per_field_co_SI==inf))=nan;
figure('units','normalized','outerposition',[0 0 1 1])

% x= solo height
for solo_param=1:5
    switch solo_param
        case 1
            x=per_field_solo_height;
            x_str='Solo field height (Hz)';
        case 2
            x=per_field_solo_nrm_height;
            x_str='Solo norm field height';
        case 3
            x=per_field_solo_size;
            x_str='Solo field size (m)';
        case 4
            x=per_field_solo_SI;
            x_str='Solo field SI (bits/spike)';
          case 5
            x=n_fields_per_cell_for_field_scatter;
            x_str='Solo n field';
            
    end
subplot(5,6,1+6*(solo_param-1))
y=per_field_co_SI;
y_str='CO per field SI (bits/spike)';
if solo_param==5
    plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
else
    plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
end

subplot(5,6,2+6*(solo_param-1))
y=per_field_co_CV;
y_str='CO per field CV';
if solo_param==5
    plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
else
    plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
end

subplot(5,6,3+6*(solo_param-1))
y=per_field_co_iqr_by_median;
y_str='CO per field iqr/median';
if solo_param==5
    plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
else
    plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
end

subplot(5,6,4+6*(solo_param-1))
y=per_field_co_modulation_depth;
y_str='CO (max-mean)/max';
if solo_param==5
    plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
else
    plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
end

subplot(5,6,5+6*(solo_param-1))
y=per_field_co_neg_modulation_depth;
y_str='CO (mean-min)/mean';
if solo_param==5
    plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
else
    plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
end

subplot(5,6,6+6*(solo_param-1))
y=per_field_co_modulation_depth_zscore_r;
y_str='zscore (max-min)/max';
if solo_param==5
    plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
else
    plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
end
end

% save figure:
file_name=fullfile(figure_folder,'scatter_per_field_box.jpg');
saveas(gcf,file_name)

%%  PLOT SUBfunctions
function plot_scatter(x,y,x_name,y_name)
[r,p]=corr(x,y,'rows','pairwise');
scatter(x,y)
xlabel(x_name)
ylabel(y_name)
title(sprintf('r=%.2f p=%.2f n=%d',r,p,length(x)))
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

legend('between cells','within a cell')
[h,p] = kstest2(x_data,shuffle);

title(sprintf('%s, p=%.2f n=%d', name,p,length(x_data)))
end

 function hist_per_field(data,x_vec,txt,xlimits,x_str)
[h,x]=hist(data,x_vec);
bar(x,h/sum(h))
n_fields=length(data);
xlabel(x_str)
ylabel('proportion')
mean_hist=nanmean(data);
median_hist=nanmedian(data);
title(sprintf('%s \n #fields=%d\nmean=%.2f meadian=%.2f',txt,n_fields,mean_hist,median_hist))
xlim(xlimits)
 end

 function plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
 scatter(x(find(id_per_field_non_signif)),y(find(id_per_field_non_signif)),2,[.5 .5 .5])
 hold on;
 scatter(x(find(id_per_field_neg_tuning)),y(find(id_per_field_neg_tuning)),5,[1 0 0],'filled')
 scatter(x(find(id_per_field_pos_tuning)),y(find(id_per_field_pos_tuning)),8,[0 1 0])
 
 % compute corrs:
 [r_non_signif,p_non_signif]=corr(x(find(id_per_field_non_signif))',y(find(id_per_field_non_signif))','rows','pairwise','Type','Spearman');
 [r_neg,p_neg]=corr(x(find(id_per_field_neg_tuning))',y(find(id_per_field_neg_tuning))','rows','pairwise','Type','Spearman');
 [r_pos,p_pos]=corr(x(find(id_per_field_pos_tuning))',y(find(id_per_field_pos_tuning))','rows','pairwise','Type','Spearman');
 %cmp_ind=find(id_per_field_pos_tuning & id_per_field_pos_tuning);
 %[r_cmp,p_cmp]=corr(x(cmp_ind),y(cmp_ind)','rows','pairwise');
 [r_all,p_all]=corr(x',y','rows','pairwise','Type','Spearman');
 title(sprintf('ALL: rho=%.2f p=%.3f',r_all,p_all))
 legend(sprintf('n.s:rho=%.2f p=%.3f',r_non_signif,p_non_signif),sprintf('neg:rho=%.2f p=%.3f',r_neg,p_neg),sprintf('pos:rho=%.2f p=%.3f',r_pos,p_pos))
 xlabel(x_str)
 ylabel(y_str)
 end
 
 function plot_scatter_per_field_with_boxplot(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)
 scatter(x(find(id_per_field_non_signif)),y(find(id_per_field_non_signif)),2,[.5 .5 .5])
 hold on;
 scatter(x(find(id_per_field_neg_tuning)),y(find(id_per_field_neg_tuning)),5,[1 0 0],'filled')
 scatter(x(find(id_per_field_pos_tuning)),y(find(id_per_field_pos_tuning)),8,[0 1 0])
 
 % compute corrs:
 [r_non_signif,p_non_signif]=corr(x(find(id_per_field_non_signif))',y(find(id_per_field_non_signif))','rows','pairwise','Type','Spearman');
 [r_neg,p_neg]=corr(x(find(id_per_field_neg_tuning))',y(find(id_per_field_neg_tuning))','rows','pairwise','Type','Spearman');
 [r_pos,p_pos]=corr(x(find(id_per_field_pos_tuning))',y(find(id_per_field_pos_tuning))','rows','pairwise','Type','Spearman');
 %cmp_ind=find(id_per_field_pos_tuning & id_per_field_pos_tuning);
 %[r_cmp,p_cmp]=corr(x(cmp_ind),y(cmp_ind)','rows','pairwise');
 [r_all,p_all]=corr(x',y','rows','pairwise','Type','Spearman');
 title(sprintf('ALL: rho=%.2f p=%.3f',r_all,p_all))
 lh=legend(sprintf('n.s:rho=%.2f p=%.3f',r_non_signif,p_non_signif),sprintf('neg:rho=%.2f p=%.3f',r_neg,p_neg),sprintf('pos:rho=%.2f p=%.3f',r_pos,p_pos))
 %set(0,'DefaultLegendAutoUpdate','off')

 xlabel(x_str)
 ylabel(y_str)
 
 % box plot:
 for cells_i=2:3
     switch cells_i
         case 1
             id_cells=id_per_field_non_signif;
             color_cells=[0.5 0.5 0.5];
         case 2
             id_cells=id_per_field_neg_tuning;
             color_cells=[1 0 0];
         case 3
             id_cells=id_per_field_pos_tuning;
             color_cells=[0 1 0];
     end
M=[y(find(id_cells))',x(find(id_cells))'];
[~,~,X] = unique(M(:,2));
C = accumarray(X,1:size(M,1),[],@(r){M(r,1)});
median_data=cellfun(@nanmedian,C,'UniformOutput',false);
extreme_prctile_data=cellfun(@(x) prctile(x,[10,90]),C,'UniformOutput',false);
prctile_data=cellfun(@(x) prctile(x,[25,75]),C,'UniformOutput',false);

plot(unique(M(:,2))+0.2*cells_i,[median_data{:}],'+','color',color_cells)
plot(repmat(unique(M(:,2)),1,2)+0.2*cells_i,reshape([prctile_data{:}],2,[])','.','color',color_cells,'markersize',3)
plot(repmat(unique(M(:,2)),1,2)+0.2*cells_i,reshape([extreme_prctile_data{:}],2,[])','.','color',color_cells,'markersize',1)
 
 end
s = get(lh, 'string'); 
set(lh, 'string', s(1:3));
 end
 function mat=convert_to_mat(cell)
    maxSize = max(cellfun(@numel,cell));
    fcn = @(x) [x nan(1,maxSize-numel(x))];
    % a. time to co
    rmat = cellfun(fcn,cell,'UniformOutput',false);
    mat = vertcat(rmat{:});
end

