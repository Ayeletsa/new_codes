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

%% initialize:
per_field_cell_dir_count=0;
fr_count_all=0;
cell_dir_per_field_count=0;

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
       
        %% ego signif cell
 
%         shuffle_struct_name = ['co_shuffling_struct_b',num2str(bat),'_d',num2str( day),'_c',num2str(cell_num),'.mat'];
%         shuffle_file_name = fullfile(co_shuffle_structs_folder,shuffle_struct_name);

        ego_signif_cells(cell_i,dir_i)=inclusion(dir_i).ego_cell;
        
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
                %2. get solo prop per field:
                per_field_solo_height(fr_count_all)=solo_field_data(field_i).peak;
                per_field_solo_nrm_height(fr_count_all)=solo_field_data(field_i).peak./nanmean(cell_co_solo_initial_analysis.solo(dir_i).PSTH_for_field_detection);
                per_field_solo_size(fr_count_all)=solo_field_data(field_i).width_href;
                per_field_solo_SI(fr_count_all)=solo_field_data(field_i).field_SI;
                 
                %3. get co prop per field:
                per_field_co_SI(fr_count_all)=per_field_data(field_i).SI;
                per_field_co_CV(fr_count_all)=per_field_data(field_i).cv;
                per_field_co_modulation_depth(fr_count_all)=per_field_data(field_i).modulation_depth;
                per_field_co_sparsity(fr_count_all)=per_field_data(field_i).sparsity;
                
                %4. get tuning width and rise time:
                if per_field_data(field_i).time_signif_field.signif_based_on_extreme_bins==1 
                    rise_and_width_data=get_per_field_tuning_width_and_rise_time(per_field_data(field_i).time_signif_field,min_dis_pos_neg);
                    id_per_field_non_signif(fr_count_all)=0;
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
                    
              
                else
                    id_per_field_non_signif(fr_count_all)=1;
                    id_per_field_pos_tuning(fr_count_all)=0;
                    id_per_field_neg_tuning(fr_count_all)=0;
                  
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
% b. compute correlations within a cell and between cells:
data=per_field_tuning_curve_relevant_cells;
[per_field_cells_corr_between,per_field_cells_corr_within,per_field_mean_corr_within]=compute_corrs_within_and_between_cells_per_field(data);


%% PLOT
%1. plot direction tuning corrs:
figure('units','normalized','outerposition',[0 0 1 1])
% ego tuning:
%-----------------
% corr between directions all cells
subplot(4,7,1:2)
plot_hists(all_cells_ego_corr_within,all_cells_ego_corr_between,'all cells egocentric corr (between directions)')

% corr between directions ego signif cells
subplot(4,7,3:4)
plot_hists(ego_signif_cells_ego_corr_within,ego_signif_cells_ego_corr_between,'ego signif cells egocentric corr (between directions)')

% Per field
%--------------------------------
% corr between per field for valid cells
subplot(4,7,8:9)
plot_hists(per_field_cells_corr_within,per_field_cells_corr_between,'per field egocentric correlations')

% scatter of isolation distance and per field mean corr:
subplot(4,7,10:11)
plot_scatter(per_field_mean_corr_within',Isolation_dis_relevant_cells','mean corr per field','Isolation distance')

%scatter of n fields and mean corr
subplot(4,7,12:13)
plot_scatter(per_field_mean_corr_within',n_fields_per_cell_relevant_cells','mean corr per field','n fields')

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
subplot(5,4,1+4*(solo_param-1))
y=per_field_co_SI;
y_str='CO per field SI (bits/spike)';
plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)

subplot(5,4,2+4*(solo_param-1))
y=per_field_co_CV;
y_str='CO per field CV';
plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)

subplot(5,4,3+4*(solo_param-1))
y=per_field_co_modulation_depth;
y_str='CO per field Modulation depth';
plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)

subplot(5,4,4+4*(solo_param-1))
y=per_field_co_sparsity;
y_str='CO per field sparsity';
plot_scatter_per_field(x,y,id_per_field_non_signif,id_per_field_pos_tuning,id_per_field_neg_tuning,x_str,y_str)

end

% save figure:
file_name=fullfile(figure_folder,'scatter_per_field.jpg');
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
 [r_non_signif,p_non_signif]=corr(x(find(id_per_field_non_signif))',y(find(id_per_field_non_signif))','rows','pairwise');
 [r_neg,p_neg]=corr(x(find(id_per_field_neg_tuning))',y(find(id_per_field_neg_tuning))','rows','pairwise');
 [r_pos,p_pos]=corr(x(find(id_per_field_pos_tuning))',y(find(id_per_field_pos_tuning))','rows','pairwise');
 %cmp_ind=find(id_per_field_pos_tuning & id_per_field_pos_tuning);
 %[r_cmp,p_cmp]=corr(x(cmp_ind),y(cmp_ind)','rows','pairwise');
 [r_all,p_all]=corr(x',y','rows','pairwise');
 title(sprintf('ALL: r=%.2f p=%.3f',r_all,p_all))
 legend(sprintf('n.s:r=%.2f p=%.3f',r_non_signif,p_non_signif),sprintf('neg:r=%.2f p=%.3f',r_neg,p_neg),sprintf('pos:r=%.2f p=%.3f',r_pos,p_pos))
 xlabel(x_str)
 ylabel(y_str)
 end
%%
%         %% OLD:
%         if  ~isempty([cell_co_solo_initial_analysis.solo(dir_i).field_height]) && e==1
%             cell_count=cell_count+1;
%             Isolation_dis(cell_count)=cell_co_solo_initial_analysis.exp_data.Isolation_dis;
%             L_Ratio(cell_count)=cell_co_solo_initial_analysis.exp_data.L_Ratio;
%             fr_count=0;
%             for fr_i=1:length([cell_co_solo_initial_analysis.solo(dir_i).field_height])
%                 if sum(isnan(cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_x_fr_per_field{1, fr_i}))<0.5*length(cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_x_fr_per_field{1, fr_i})
%                     fr_count=fr_count+1;
%                     fr_count_all=fr_count_all+1;
%                     per_field_tuning_curve_all_cells{cell_count,fr_count}=cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_x_fr_per_field{1, fr_i}';
%                     firing_rate=cell_co_solo_initial_analysis.co(dir_i).firing_rate;
%                     solo=cell_co_solo_initial_analysis.solo(dir_i);
%                     per_field_SI(fr_count_all)=firing_rate.information_per_spike_per_field{fr_i};
%                     per_field_CV(fr_count_all)=firing_rate.cv{fr_i};
%                     per_field_modulation_depth(fr_count_all)=firing_rate.modulation_depth{fr_i};
%                     per_field_sparsity(fr_count_all)=firing_rate.sparsity{fr_i};
%                     per_field_n_spikes(fr_count_all)=firing_rate.number_of_spikes_per_field{fr_i};
%                     per_field_height(fr_count_all)=solo.field_height(fr_i);
%                     per_field_nrm_height(fr_count_all)=solo.field_height_norm_by_mean(fr_i);
%                     per_field_size(fr_count_all)=solo.field_size(fr_i);
%                     
%                     
%                     
%                     
%                     per_field_signif(cell_count)=0;% it will change if it will be true
%                     for fr_i=1:length(firing_rate.time_signif_field_new_smooth)
%                         
%                         if iscell(firing_rate.time_signif_field_new_smooth{fr_i}) & firing_rate.number_of_spikes_per_field{fr_i}>min_n_spike_per_field
%                             
%                             for smooth_i=1:n_smooth
%                                 r=firing_rate.time_fr_per_field_new_smooth_smooth_window{fr_i}(smooth_i,:);
%                                 non_nan_ind_r=find(~isnan(r));
%                                 [ind_length,~,~]=find_length_of_consecutive_ind(non_nan_ind_r,length(r));
%                                 
%                                 if firing_rate.time_signif_field_new_smooth{fr_i}.signif_based_on_extreme_bins==1 & ind_length>min_r_length_per_field*length(r)
%                                     
%                                     per_field_signif(cell_count)=1;
%                                     neg_width=firing_rate.time_signif_field_new_smooth{fr_i}.width_neg_interp;
%                                     neg_width(isnan(neg_width))=[];
%                                     pos_width=firing_rate.time_signif_field_new_smooth{fr_i}.width_pos_interp;
%                                     pos_width(isnan(pos_width))=[];
%                                     pos_rise_time=firing_rate.time_signif_field_new_smooth{fr_i}.pos_rise_time_interp;
%                                     pos_rise_time(isnan(pos_rise_time))=[];
%                                     neg_rise_time=firing_rate.time_signif_field_new_smooth{fr_i}.neg_rise_time_interp;
%                                     neg_rise_time(isnan(neg_rise_time))=[];
%                                     
%                                     comb_width=[];
%                                     pos_comb=[];
%                                     neg_comb=[];
%                                     
%                                     if ~isempty(pos_width) & ~isempty(neg_width) % if there are both positive and negative
%                                         x_pos=firing_rate.time_signif_field_new_smooth{fr_i}.width_line_x_pos_interp;
%                                         x_neg=firing_rate.time_signif_field_new_smooth{fr_i}.width_line_x_neg_interp;
%                                         x_pos(isnan(x_pos))=[];x_neg(isnan(x_neg))=[];
%                                         % check if there are clsoe to each
%                                         % other:
%                                         if abs(min(x_pos(:))-max(x_neg(:)))<=min_dis_pos_neg | abs(max(x_pos(:))-min(x_neg(:)))<=min_dis_pos_neg
%                                             %TO DO: to correct it to work
%                                             %properly on all type of tuning
%                                             %curve if we will use it!
%                                             comb_width=[max([x_pos(:);x_neg(:)])-min([x_pos(:);x_neg(:)])];
%                                             pos_comb=pos_width;
%                                             neg_comb=neg_width;
%                                             pos_width=[];
%                                             neg_width=[];
%                                             neg_rise_time=[];
%                                             pos_rise_time=[];
%                                         end
%                                     end
%                                     all_pos_width=[all_pos_width,pos_width];
%                                     all_neg_width=[all_neg_width,neg_width];
%                                     all_compound_width=[all_compound_width,comb_width];
%                                     all_compound_pos_width=[all_compound_pos_width,pos_comb];
%                                     all_compound_neg_width=[all_compound_neg_width,neg_comb];
%                                     all_pos_rise_time=[all_pos_rise_time,pos_rise_time];
%                                     all_neg_rise_time=[all_neg_rise_time,neg_rise_time];
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %%
% % take only cells with more than one field to calculat:
% [r c]=find(~cellfun(@isempty,per_field_tuning_curve_all_cells));
% ind=zeros(size(per_field_tuning_curve_all_cells));
% ind(sub2ind(size(per_field_tuning_curve_all_cells),r,c))=1;
% cells_to_remove=sum(ind,2)<=1;
% L_Ratio_relevant_cells=L_Ratio(~cells_to_remove);
% Isolation_dis_relevant_cells=Isolation_dis(~cells_to_remove);
% per_field_tuning_curve_relevant_cells=per_field_tuning_curve_all_cells(~cells_to_remove,:);
% all_cells_corr=[];
% other_cells_corr=[];
% for cell_i=1:length(per_field_tuning_curve_relevant_cells)
%     cell_vec=zeros(length(per_field_tuning_curve_relevant_cells),1);
%     cell_vec(cell_i)=1;
%     mat_1=[per_field_tuning_curve_relevant_cells{cell_i,:}];
%     cell_corr_values=corr(mat_1,mat_1,'rows','pairwise');
%     cell_corr_values=triu(cell_corr_values,1);
%     cell_corr_values=cell_corr_values(find(cell_corr_values~=0));
%     all_cells_corr=[all_cells_corr; cell_corr_values];
%     mean_corr_cell_value(cell_i)=nanmean(cell_corr_values);
%     
%     other_cells_mat=[per_field_tuning_curve_relevant_cells{~cell_vec,:}];
%     
%     other_cell_corr_values=corr(mat_1,other_cells_mat,'rows','pairwise');
%     other_cells_corr=[other_cells_corr;other_cell_corr_values(:)];
% end
% 
% % only signif per field:
% relevant_cells=~cells_to_remove & per_field_signif';
% per_field_tuning_curve_relevant_cells=per_field_tuning_curve_all_cells(relevant_cells,:);
% cells_with_signif_bin_corr=[];
% other_cells_with_signif_bin_cells_corr=[];
% for cell_i=1:length(per_field_tuning_curve_relevant_cells)
%     cell_vec=zeros(length(per_field_tuning_curve_relevant_cells),1);
%     cell_vec(cell_i)=1;
%     mat_1=[per_field_tuning_curve_relevant_cells{cell_i,:}];
%     cell_corr_values=corr(mat_1,mat_1,'rows','pairwise');
%     cell_corr_values=triu(cell_corr_values,1);
%     cells_with_signif_bin_corr=[cells_with_signif_bin_corr; cell_corr_values(find(cell_corr_values~=0))];
%     
%     other_cells_mat=[per_field_tuning_curve_relevant_cells{~cell_vec,:}];
%     
%     other_cell_corr_values=corr(mat_1,other_cells_mat,'rows','pairwise');
%     other_cells_with_signif_bin_cells_corr=[other_cells_with_signif_bin_cells_corr;other_cell_corr_values(:)];
% end
% 
% % egocentric signif cells tuning
% signif_ego_corr_cell=[];
% signif_ego_corr_other_cells=[];
% any_dir_signif_cells=sum(ego_signif_cells,2);
% ego_signif_tuning=all_ego_tuning(any_dir_signif_cells>0,:,:);
% for cell_i=1:size(ego_signif_tuning,1)
%     
%     mat_1=squeeze(ego_signif_tuning(cell_i,:,:));
%     cell_corr_values=corr(mat_1(1,:)',mat_1(2,:)','rows','pairwise');
%     signif_ego_corr_cell=[signif_ego_corr_cell,cell_corr_values];
%     
%     cell_vec=zeros(size(ego_signif_tuning,1),1);
%     cell_vec(cell_i)=1;
%     other_cells=ego_signif_tuning(~cell_vec,:,:);
%     other_cells=[squeeze(other_cells(:,1,:));squeeze(other_cells(:,2,:))];
%     other_cell_corr_values=corr(mat_1(1,:)',other_cells','rows','pairwise');
%     signif_ego_corr_other_cells=[signif_ego_corr_other_cells;other_cell_corr_values(:)];
%     other_cell_corr_values=corr(mat_1(2,:)',other_cells','rows','pairwise');
%     signif_ego_corr_other_cells=[signif_ego_corr_other_cells;other_cell_corr_values(:)];
% end
% 
% % all egocentric tuning
% all_ego_corr_cell=[];
% all_ego_corr_other_cells=[];
% for cell_i=1:length(all_ego_tuning)
%     mat_1=squeeze(all_ego_tuning(cell_i,:,:));
%     cell_corr_values=corr(mat_1(1,:)',mat_1(2,:)','rows','pairwise');
%     all_ego_corr_cell=[all_ego_corr_cell,cell_corr_values];
%     
%     cell_vec=zeros(length(all_ego_tuning),1);
%     cell_vec(cell_i)=1;
%     other_cells=all_ego_tuning(~cell_vec,:,:);
%     other_cells=[squeeze(other_cells(:,1,:));squeeze(other_cells(:,2,:))];
%     other_cell_corr_values=corr(mat_1(1,:)',other_cells','rows','pairwise');
%     all_ego_corr_other_cells=[all_ego_corr_other_cells;other_cell_corr_values(:)];
%     other_cell_corr_values=corr(mat_1(2,:)',other_cells','rows','pairwise');
%     all_ego_corr_other_cells=[all_ego_corr_other_cells;other_cell_corr_values(:)];
% end
% %% PLOT
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% subplot(3,2,1)
% plot_hists(all_ego_corr_cell,all_ego_corr_other_cells,'all cells egocentric corr (between directions)')
% 
% subplot(3,2,2)
% plot_hists(signif_ego_corr_cell,signif_ego_corr_other_cells,'ego signif cells egocentric corr (between directions)')
% 
% subplot(3,2,3)
% plot_hists(all_cells_corr,other_cells_corr,'per field egocentric correlations')
% 
% subplot(3,2,4)
% plot_hists(cells_with_signif_bin_corr,other_cells_with_signif_bin_cells_corr,'per field egocentric correlations- cells with signif bin in at least one field')
% 
% subplot(3,2,5)
% plot_scatter(mean_corr_cell_value',L_Ratio_relevant_cells','mean corr per field','L Ratio')
% 
% subplot(3,2,6)
% plot_scatter(mean_corr_cell_value',Isolation_dis_relevant_cells','mean corr per field','Isolation distance')
% 
% % save figure:
% file_name=fullfile(figure_folder,'hist_corr_directions_and_per_field.jpg');
% saveas(gcf,file_name)
% %% Figure of per field tuning with
% figure('units','normalized','outerposition',[0 0 1 1])
% smooth_vec=[3,5,7];
% 
% n_plots=7;
% for smooth_i=1:n_smooth
%     x_vec=(0:0.2:4);
%     xlimits=[0 4];
%     smooth_wind=smooth_vec(smooth_i);
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+1)
%     data=all_pos_width./1e6;
%     txt='Positive fields:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+2)
%     data=all_neg_width./1e6;
%     txt='Negative fields:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+3)
%     data=all_compound_width./1e6;
%     txt='Compound fields:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+4)
%     data=all_compound_pos_width./1e6;
%     txt='Pos from conpound:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+5)
%     data=all_compound_neg_width./1e6;
%     txt='Neg from compound:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
%     x_vec=(0:0.1:2);
%     xlimits=[0 2];
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+6)
%     data=all_pos_rise_time./1e6;
%     txt='Pos rise time:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
%     subplot(n_smooth,n_plots,n_plots*(smooth_i-1)+7)
%     data=all_neg_rise_time./1e6;
%     txt='Neg rise time:';
%     hist_per_field(data,x_vec,txt,smooth_wind,xlimits)
%     
% end
% 
% file_name=fullfile(figure_folder,'hist_per_field_tuning_width.jpg');
% saveas(gcf,file_name)
% %%
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(4,4,1)
% plot_scatter(per_field_n_spikes',per_field_SI','# spikes','SI')
% subplot(4,4,2)
% plot_scatter(per_field_n_spikes',per_field_CV','# spikes','CV')
% subplot(4,4,3)
% plot_scatter(per_field_n_spikes',per_field_modulation_depth','# spikes','modulation depth')
% subplot(4,4,4)
% plot_scatter(per_field_n_spikes',per_field_sparsity','# spikes','sparsity')
% 
% subplot(4,4,5)
% plot_scatter(per_field_height',per_field_SI','field height','SI')
% subplot(4,4,6)
% plot_scatter(per_field_height',per_field_CV','field height','CV')
% subplot(4,4,7)
% plot_scatter(per_field_height',per_field_modulation_depth','field height','modulation depth')
% subplot(4,4,8)
% plot_scatter(per_field_height',per_field_sparsity','field height','sparsity')
% 
% subplot(4,4,9)
% plot_scatter(per_field_nrm_height',per_field_SI','norm field height','SI')
% subplot(4,4,10)
% plot_scatter(per_field_nrm_height',per_field_CV','norm field height','CV')
% subplot(4,4,11)
% plot_scatter(per_field_nrm_height',per_field_modulation_depth','norm field height','modulation depth')
% subplot(4,4,12)
% plot_scatter(per_field_nrm_height',per_field_sparsity','norm field height','sparsity')
% 
% 
% subplot(4,4,13)
% plot_scatter(per_field_size',per_field_SI','field size','SI')
% subplot(4,4,14)
% plot_scatter(per_field_size',per_field_CV','field size','CV')
% subplot(4,4,15)
% plot_scatter(per_field_size',per_field_modulation_depth','field size','modulation depth')
% subplot(4,4,16)
% plot_scatter(per_field_size',per_field_sparsity','field size','sparsity')
% 
% % save figure:
% file_name=fullfile(figure_folder,'per_field_properties_corr.jpg');
% saveas(gcf,file_name)
% 
% 
