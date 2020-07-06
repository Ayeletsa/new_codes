%% session b analysis:
clear
% params
min_co_dir=10;
min_n_co_spikes=30;
load('D:\Ayelet\2bat_proj\Analysis\new_code\params\per_field_params.mat')
%% arrange data:
inclusion_cells_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\inclusion_cells_struct';

relevant_days={'20191210','20191212','20191213'};

dir_data_session_a='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_data_session_b='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\session_b\';
dir_info=dir(dir_data_session_a);
dir_info_names_session_a={dir_info.name};
dir_info_names_session_a([dir_info.isdir])=[];
dir_info=dir(dir_data_session_b);
dir_info_names_session_b={dir_info.name};
dir_info_names_session_b([dir_info.isdir])=[];

% 1. find days ind:
%-------------------------
day_ind=regexp(dir_info_names_session_a{1},'day_');
day_names=cellfun(@(a) a(day_ind+4:day_ind+11),dir_info_names_session_a,'UniformOutput',false);
Match=cellfun(@(a) find(contains(a,relevant_days)) , day_names, 'UniformOutput', 0);
relevant_day_ind=find(cellfun(@(c) ~isempty(c),Match));

%2. find cell num in days:
%-------------------------
cell_ind=regexp(dir_info_names_session_a{1},'cell_');
cell_nums=cellfun(@(a) a(cell_ind+5:cell_ind+7),dir_info_names_session_a,'UniformOutput',false);
cell_nums=cellfun(@(c) str2num(c),cell_nums,'UniformOutput',false);
cell_nums_in_day=cell_nums(relevant_day_ind);
cell_nums_in_day=cell2mat(cell_nums_in_day);
% % find cells that are signif ego in session a
% %-------------------------
% dir_info=dir(inclusion_cells_folder_name);
% dir_info_names={dir_info.name};
% dir_info_names([dir_info.isdir])=[];
% dir_info_names=dir_info_names(2:end);
% %cell_ind=regexp(dir_info_names{2},'cell_');
% cell_ind=16;
% cell_nums_inclusion=cellfun(@(a) a(cell_ind:cell_ind+2),dir_info_names,'UniformOutput',false);
% cell_nums_inclusion=cellfun(@(c) str2num(c),cell_nums_inclusion,'UniformOutput',false);
% [~, r]=ismember(cell2mat(cell_nums_in_day),cell2mat(cell_nums_inclusion))
% %Match=cellfun(@(a) find(contains(a,strsplit(num2str(cell_nums_in_day)))) , dir_info_names, 'UniformOutput', 0);
% %r=find(cellfun(@(c) ~isempty(c),Match));
% 
% %loda data
% all_cells_data=struct();
% for cell_i=1:length(r)
%     cell_file_name=dir_info_names{r(cell_i)};
%     cell_ind=regexp(cell_file_name,'cell_');
%     
%     file_names=fullfile(inclusion_cells_folder_name,cell_file_name);
%     cell_num=str2num(cell_file_name(cell_ind+5:cell_ind+7));
%     load(file_names)
%     all_cells_data(cell_i).('ego_cell') = [inclusion.('ego_cell')];
%     all_cells_data(cell_i).('place_cell') = [inclusion.('place_cell')];
%     all_cells_data(cell_i).('cell_num') =cell_num;
% end

% run over place cells:
% place_cells_ind=find(sum(reshape([all_cells_data.('place_cell')],2,[])));
% place_cells_num=[all_cells_data(place_cells_ind).cell_num];
% place=1;
% [solo_fr_place_cells,co_fr_place_cells,allo_ego_map_place_cells]=get_cell_sessions_data(place_cells_num,dir_info_names_session_a,dir_data_session_a,dir_info_names_session_b,dir_data_session_b,inclusion_cells_folder_name,min_co_dir,min_n_co_spikes);


% run over ego cells:
% ego_cells_ind=find(sum(reshape([all_cells_data.('ego_cell')],2,[])));
% ego_cells_num=[all_cells_data(ego_cells_ind).cell_num];
% place=0;
[solo_fr_ego_cells,co_fr_ego_cells,allo_ego_map_ego_cells,all_per_field_tuning]=get_cell_sessions_data(cell_nums_in_day,dir_info_names_session_a,dir_data_session_a,dir_info_names_session_b,dir_data_session_b,inclusion_cells_folder_name,min_co_dir,min_n_co_spikes,min_r_length_per_field,min_n_spike_per_field);

%% plot hist of corrs for allo and ego tuning
figure('units','normalized','outerposition',[0 0 1 1])

% Ego cells:
%--------------------------
% ego tuning
cells=1;
[corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(co_fr_ego_cells);
subplot(2,5,1)
plot_corr=1;
name='corr ego tuning for ego cells';
x_data=corr_cell;
shuffle=corr_all_other_cells;
plot_hists(x_data,shuffle,name,plot_corr,cells)

subplot(2,5,2)
plot_corr=0;
name='delta peak ego tuning for ego cells';
x_data=delta_peak_cell;
shuffle=delta_peak_all_other_cells;
plot_hists(x_data,shuffle,name,plot_corr,cells)

%place tuning
[corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(solo_fr_ego_cells);
subplot(2,5,3)
plot_corr=1;
name='corr allo tuning for ego cells';
x_data=corr_cell;
shuffle=corr_all_other_cells;
plot_hists(x_data,shuffle,name,plot_corr,cells)

subplot(2,5,4)
plot_corr=0;
name='delta peak allo tuning for ego cells';
x_data=delta_peak_cell;
shuffle=delta_peak_all_other_cells;
plot_hists(x_data,shuffle,name,plot_corr,cells)

%2D tuning
[corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(allo_ego_map_ego_cells);
plot_corr=1;
subplot(2,5,5)
name='corr 2D tuning for ego cells';
x_data=corr_cell;
shuffle=corr_all_other_cells;
plot_hists(x_data,shuffle,name,plot_corr,cells)

%2D tuning
[corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(all_per_field_tuning);
plot_corr=1;
subplot(2,5,6)
name='corr per field ego cells';
x_data=corr_cell;
shuffle=corr_all_other_cells;
cells=0;
plot_hists(x_data,shuffle,name,plot_corr,cells)

% % Place cells
% %------------------
% % ego tuning
% plot_corr=1;
% [corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(co_fr_place_cells);
% subplot(2,5,6)
% name='corr ego tuning for place cells';
% x_data=corr_cell;
% shuffle=corr_all_other_cells;
% plot_hists(x_data,shuffle,name,plot_corr)
% 
% plot_corr=0;
% subplot(2,5,7)
% name='delta peak ego tuning for place cells';
% x_data=delta_peak_cell;
% shuffle=delta_peak_all_other_cells;
% plot_hists(x_data,shuffle,name,plot_corr)
% 
% 
% %place tuning
% [corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(solo_fr_place_cells);
% plot_corr=1;
% subplot(2,5,8)
% name='corr allo tuning for place cells';
% x_data=corr_cell;
% shuffle=corr_all_other_cells;
% plot_hists(x_data,shuffle,name,plot_corr)
% 
% plot_corr=0;
% subplot(2,5,9)
% name='delta peak allo tuning for place cells';
% x_data=delta_peak_cell;
% shuffle=delta_peak_all_other_cells;
% plot_hists(x_data,shuffle,name,plot_corr)
% 
% 
% %2D tuning
% [corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(allo_ego_map_place_cells);
% 
% subplot(2,5,10)
% name='corr 2D tuning for place cells';
% x_data=corr_cell;
% shuffle=corr_all_other_cells;
% plot_corr=1;
% plot_hists(x_data,shuffle,name,plot_corr)
% 
%% save
direcoty='D:\Ayelet\2bat_proj\Analysis\new_code\figures\population\';
file_name=fullfile(direcoty,'sessiona_b_corrs.jpg');
saveas(gcf,file_name);

%% functions
function [solo_fr,co_fr,allo_ego_map,all_per_field_tuning]=get_cell_sessions_data(cell_num,dir_info_names_session_a,dir_data_session_a,dir_info_names_session_b,dir_data_session_b,inclusion_cells_folder_name,min_co_dir,min_n_co_spikes,min_r_length_per_field,min_n_spike_per_field)
count=0;
cell_dir_count=0;
for cell_i=cell_num
    %load both sessions data:
    cell_ind=regexp(dir_info_names_session_a{1},'cell_');
    cell_nums=cellfun(@(a) a(cell_ind+5:cell_ind+7),dir_info_names_session_a,'UniformOutput',false);
    Match=cellfun(@(a) find(contains(a,strsplit(num2str(cell_i)))) , cell_nums, 'UniformOutput', 0);
    r=find(cellfun(@(c) ~isempty(c),Match));
    cell_data_session_a=load(fullfile(dir_data_session_a,dir_info_names_session_a{r}));
    
    cell_ind=regexp(dir_info_names_session_b{1},'cell_');
    cell_nums=cellfun(@(a) a(cell_ind+5:cell_ind+7),dir_info_names_session_b,'UniformOutput',false);
    Match=cellfun(@(a) find(contains(a,strsplit(num2str(cell_i)))) , cell_nums, 'UniformOutput', 0);
    r=find(cellfun(@(c) ~isempty(c),Match));
    if ~isempty(r)
        cell_data_session_b=load(fullfile(dir_data_session_b,dir_info_names_session_b{r}));
        
        file_name=fullfile(inclusion_cells_folder_name,['inclusion_cell_',num2str(cell_i),'.mat']);
        load(file_name)
       % if place==1
       relevant_session=[inclusion.ego_cell]; 
       for dir_i=1:2
           per_field=cell_data_session_a.cell_co_solo_initial_analysis.co(dir_i).per_field_href;
           if ~isempty(fieldnames(per_field)) & inclusion(dir_i).place_cell
               for field=1:length(per_field)
                   if ~isempty(per_field(field).dis_signif_field)
                       if per_field(field).dis_signif_field.signif_based_on_extreme_bins
                           relevant_session(dir_i)=1;
                       end
                   end
               end
           end
       end
       

         for dir_i=1:2
             cond(1)=cell_data_session_b.cell_co_solo_initial_analysis.exp_data.num_co_per_dir(dir_i)>min_co_dir;
             cond(2)=cell_data_session_b.cell_co_solo_initial_analysis.co(1).info.n_spikes  >min_n_co_spikes;
             valid=all(cond);
            if relevant_session(dir_i)==1 & valid
              
                cell_dir_count=cell_dir_count+1;
                solo_fr(:,1,cell_dir_count)=cell_data_session_a.cell_co_solo_initial_analysis.solo(dir_i).x_pos_firing_rate{1, 1}' ;
                solo_fr(:,2,cell_dir_count)=cell_data_session_b.cell_co_solo_initial_analysis.solo(dir_i).x_pos_firing_rate{1, 1}' ;
                co_fr(:,1,cell_dir_count)=cell_data_session_a.cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:)' ;
                co_fr(:,2,cell_dir_count)=cell_data_session_b.cell_co_solo_initial_analysis.co(dir_i).firing_rate.dis_m(1,:)' ;
                allo_ego_map(:,1,cell_dir_count)=cell_data_session_a.cell_co_solo_initial_analysis.co(dir_i).firing_rate.field_density_smoothed_XY_with_NaN(:);
                allo_ego_map(:,2,cell_dir_count)=cell_data_session_b.cell_co_solo_initial_analysis.co(dir_i).firing_rate.field_density_smoothed_XY_with_NaN(:);
                
                % per field:
                 solo_data_session_a=cell_data_session_a.cell_co_solo_initial_analysis.solo;

                fields=solo_data_session_a(dir_i).fields;
                if ~isempty(fields)
                    field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2);
                    spikes=cell_data_session_b.cell_co_solo_initial_analysis.co(dir_i).spikes;
                    bsp=cell_data_session_b.cell_co_solo_initial_analysis.co(dir_i).bsp;
                    solo_data=cell_data_session_b.cell_co_solo_initial_analysis.solo;
                    per_field_params_file_name='D:\Ayelet\2bat_proj\Analysis\new_code\params\per_field_params.mat';
                    [per_field_href,per_field_tunings_corrs]=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(dir_i));
                    
                end
                  
%                     field_edges=reshape([fields.edges_href],2,length([fields.edges_href])/2);
%                     
%                     per_field_params_file_name
%                     [per_field_href,per_field_tunings_corrs]=per_field_analysis(field_edges,spikes,bsp,per_field_params_file_name,solo_data(ii_dir));
%                     
                fields_session_a=cell_data_session_a.cell_co_solo_initial_analysis.solo(dir_i).fields;
                %fields_session_b=cell_data_session_b.cell_co_solo_initial_analysis.solo(dir_i).fields;
                for field_i=1:length(fields_session_a)
                    %sessio_b_field_locs=[fields_session_b.loc];
                    %field_edged_session_a=reshape([fields_session_a.edges_href],2,length([fields_session_a.edges_href])/2);
                    %session_b_field_ind=find(sessio_b_field_locs>field_edged_session_a(1,field_i) & sessio_b_field_locs<field_edged_session_a(2,field_i));
                    %if ~isempty(session_b_field_ind)
                    tuning_session_a=cell_data_session_a.cell_co_solo_initial_analysis.co(dir_i).per_field_href(field_i).tuning_dis_x_fr_per_field  ;
                    %find if per field is valid
                    r=tuning_session_a;
                    [ind_length,~,~]=find_length_of_consecutive_ind(find(~isnan(r)),length(r));
                    a=max(ind_length)>=min_r_length_per_field*length(r); %enough data that is not nan in tuning curve
                    b=cell_data_session_a.cell_co_solo_initial_analysis.co(dir_i).per_field_href(field_i).number_of_spikes_per_field>=min_n_spike_per_field;% enough spikes
                    per_field_cond=a & b;
                    if per_field_cond
                    tuning_session_b=per_field_href(field_i).tuning_dis_x_fr_per_field  ;
                    count=count+1;
                    all_per_field_tuning(:,1,count)=tuning_session_a;
                    all_per_field_tuning(:,2,count)=tuning_session_b;

                    end
                end
            end
            
        end
    end
end

end
function [corr_cell, corr_all_other_cells,delta_peak_cell,delta_peak_all_other_cells,delta_trough_cell,delta_trough_all_other_cells]=compute_corr_of_all_other_cells(firing_rate)
% compute_corrs:
corr_all_other_cells=corr(squeeze(firing_rate(:,1,:)),squeeze(firing_rate(:,2,:)),'rows','pairwise');
idx = logical(eye(size(corr_all_other_cells)));
corr_cell=corr_all_other_cells(idx);
other_cell_ind=triu(~idx);
corr_all_other_cells=corr_all_other_cells(other_cell_ind);

%compute delta peak
[val, peaks_loc]=nanmax(firing_rate);
peaks_loc=squeeze(peaks_loc);
val=squeeze(val);
peaks_loc(:,[find(sum(isnan(val)))])=[];
all_delta=abs(peaks_loc(1,:)'-repmat(peaks_loc(2,:),length(peaks_loc),1));
idx = logical(eye(size(all_delta)));
delta_peak_cell=all_delta(idx);
delta_peak_all_other_cells=all_delta(~idx);

%compute delta peak
[val, peaks_loc]=nanmin(firing_rate);
peaks_loc=squeeze(peaks_loc);
val=squeeze(val);
peaks_loc(:,[find(sum(isnan(val)))])=[];
all_delta=abs(peaks_loc(1,:)'-repmat(peaks_loc(2,:),length(peaks_loc),1));
idx = logical(eye(size(all_delta)));
delta_trough_cell=all_delta(idx);
delta_trough_all_other_cells=all_delta(~idx);

end

function plot_hists(x_data,shuffle,name,plot_corr,cells)
if plot_corr==1
    x_bin=-1:0.1:1;
    [h,x]=hist(shuffle,x_bin);
    b1=bar(x,h/sum(h));
    b1.FaceAlpha = 0.5;
    hold on
    [h,x]=hist(x_data,x_bin);
    b2=bar(x,h/sum(h));
    b2.FaceAlpha = 0.5;
else
    [h,x]=hist(shuffle);
    b1=bar(x,h/sum(h));
    b1.FaceAlpha = 0.5;
    hold on
    [h,x]=hist(x_data);
    b2=bar(x,h/sum(h));
    b2.FaceAlpha = 0.5;
end
if  plot_corr==1
    xlabel('r (session a,session b)')
else
    xlabel('delta bin')
end
ylabel('proportion')

legend('between cells','within a cell')
[h,p_ks] = kstest2(x_data,shuffle);
[p_wx,h] = ranksum(x_data,shuffle);
if cells
title(sprintf('%s,\n p(ks)=%.2f p(w)=%.2f\n #n cells*dir=%d', name,p_ks,p_wx,length(find(~isnan(x_data)))))
else
    title(sprintf('%s,\n p(ks)=%.2f p(w)=%.2f\n #n fields=%d', name,p_ks,p_wx,length(find(~isnan(x_data)))))

end
end
