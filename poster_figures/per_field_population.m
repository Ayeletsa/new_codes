close all
clear
clc

cells_to_pick=[5,9:11,13,14,20,21,23,24,30,31,37,39:42,44,46,...
    49:51,53:55,61,69:71,74:75,81:86,...
    213,215,216,217,218,225,226,238,240,241,243,244,252,253,257,261,...
    262,266,267,270,272,273,275,276,277,278,282,284,290,291,292,294,295,296,297,298];
fig_dir_out='D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\';
behavior_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
files = dir(behavior_structs_folder);
behavior_struct_names = {files.name};
co_shuffle_structs_folder = 'D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_shuffling_struct';
files = dir(co_shuffle_structs_folder);
co_shuffle_struct_names = {files.name};

%% parameters:
%modulation_depth_criterion=[0.5, 0.25]; %take only cells with modulation depth higher than parameter
z_score_lim=9;
CV_threshold=[0.4];
bins_to_remove=[1]; %remove bins from the sides for the ordering acording to pick/trough
norm_acording_to_other_field=[0];
thresh_number_of_nan=0.5;
%initialize:
centers=cell(2,1);
ignore_extreme_bins=1;
c_max=[];
%% figure prop


figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [24 21]; % ~A4 size
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');
line_width=3;
ego_color=[1  0 1];



% z- score mat:
ax_plot{1}=axes('position',[0.08 0.6 0.2 0.2]);%zscore mat
%trough
ax_plot{2}=axes('position',[0.38 0.6 0.2 0.2]);%zscore mat


%significant bins:
ax_plot{3}=axes('position',[0.08 0.15 0.2 0.2]);%mat
ax_plot{4}=axes('position',[0.08 0.35 0.2 0.1]);%hist
%trough
ax_plot{5}=axes('position',[0.38 0.15 0.2 0.2]);%mat
ax_plot{6}=axes('position',[0.38 0.35 0.2 0.1]);%hist
%% create mat
count=0;
pos_signif={};
neg_signif={};
per_field_tuning_dir=cell(2,1);
per_field_z_score=cell(2,1);
for ii_cell=cells_to_pick
    
        %% load data
       cell_num=ii_cell;
    str=sprintf('cell_%d.mat',cell_num);
    struct_ind =find(~cellfun(@isempty,(strfind(behavior_struct_names,str))));
    struct_name=behavior_struct_names{struct_ind};
    file_name = fullfile(behavior_structs_folder,struct_name);
    load(file_name);
    if cell_num>86
    behavior_struct=cell_co_solo_initial_analysis;
    else
            load(file_name);

    end
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    cell_num=behavior_struct.exp_data.cell_num;
    number_of_dis_bins=length(behavior_struct.co(1).firing_rate.dis_X_bins_vector_of_centers{1,1});
    
    for ii_dir=1:2
        if length(behavior_struct.solo(ii_dir).field_center)>0
            all_per_filed=cell2mat({behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1, :}});
            all_per_filed=reshape(all_per_filed,number_of_dis_bins,[])';
            
            
            
            if  norm_acording_to_other_field==1
                %remove min
                min_per_field=min(all_per_filed(:));
                all_per_filed=all_per_filed-min_per_field;
                %norm to max
                max_per_filed=max(all_per_filed(:));
                all_per_filed=all_per_filed./max_per_filed;
                
            end
            
            for field_i=1:length(behavior_struct.solo(ii_dir).field_center)
                r=all_per_filed(field_i,:);
                if isnan(behavior_struct.co(ii_dir).firing_rate.information_per_spike_per_field{1, field_i})
                    continue
                else
                    cv_value=behavior_struct.co(ii_dir).firing_rate.signif_field{1, field_i}.cv;
                end
                % criterion for fields:
                % 1. if this is more than 50% nan ignore
                if sum(isnan(r))>=thresh_number_of_nan*number_of_dis_bins
                    continue
                    
                    %2. CV
                    
                    %  elseif cv<CV_threshold(CV_threshold_i)
                elseif cv_value<CV_threshold | behavior_struct.co(ii_dir).firing_rate.signif_field{field_i}.signif_based_on_CV~=1
                    
                    continue
                else
                    centers{ii_dir}=[centers{ii_dir}, behavior_struct.solo(ii_dir).field_center(field_i)];
                    per_field_tuning_dir{ii_dir}(end+1,:)=r;
                    per_field_z_score{ii_dir}(end+1,:)=behavior_struct.co(ii_dir).firing_rate.signif_field{1, field_i}.zscore_data ;
                    
                    
                    
                    %for signif bins:
                    if behavior_struct.co(ii_dir).firing_rate.signif_field{1,field_i}.signif_based_on_extreme_bins==1
                        count=count+1;
                        pos_signif{count}= behavior_struct.co(ii_dir).firing_rate.signif_field{1,field_i}.pos_signif;
                        neg_signif{count}= behavior_struct.co(ii_dir).firing_rate.signif_field{1,field_i}.neg_signif;
                        r=behavior_struct.co(ii_dir).firing_rate.dis_x_fr_per_field{1,field_i};
                        [~, max_pos{count}]=max(r);
                        [~, min_pos{count}]=min(r);
                    end
                end
                
            end
        end
        
        
        
    end
    
end
dis_bins=behavior_struct.co(1).firing_rate.dis_X_bins_vector_of_centers{1,1};
%positions of peak and min
max_pos=cell2mat(max_pos);
min_pos=cell2mat(min_pos);
%sort by max/min
[~, max_pos_sorted_ind]=sort(max_pos);
[~, min_pos_sorted_ind]=sort(min_pos);



%% significant bins:
count=0;
axes(ax_plot{3})
hold on;
for ii_field=max_pos_sorted_ind
    count=count+1;
    pos_signif_cell=pos_signif{ii_field};
    neg_signif_cell=neg_signif{ii_field};
    if ~isempty(pos_signif_cell)
        plot(dis_bins(pos_signif_cell),count*ones(1,length(pos_signif_cell)),'g.')
    end
    if ~isempty(neg_signif_cell)
        plot(dis_bins(neg_signif_cell),count*ones(1,length(neg_signif_cell)),'r.')
    end
end
%title('Sorted by peak', 'FontSize', 11, 'FontWeight', 'Bold');
xlabel('Distance between bats (m)')
ylabel('fields')
ylim([0 count+1])


count=0;
axes(ax_plot{5})
hold on;
for ii_field=min_pos_sorted_ind
    count=count+1;
    pos_signif_cell=pos_signif{ii_field};
    neg_signif_cell=neg_signif{ii_field};
    if ~isempty(pos_signif_cell)
        plot(dis_bins(pos_signif_cell),count*ones(1,length(pos_signif_cell)),'g.')
    end
    if ~isempty(neg_signif_cell)
        plot(dis_bins(neg_signif_cell),count*ones(1,length(neg_signif_cell)),'r.')
    end
end
%title('Sorted by trough', 'FontSize', 11, 'FontWeight', 'Bold');
xlabel('Distance between bats (m)')
ylabel('fields')
ylim([0 count+1])

%% hist of significant bins
pos_bins=dis_bins(cell2mat(pos_signif));
neg_bins=dis_bins(cell2mat(neg_signif));

axes(ax_plot{4})
h=hist(pos_bins,dis_bins);
bar(dis_bins,(h./sum(h))*100,'g')
title(sprintf('Distribution of significant bins\n that are higher than shuffle'))
ylabel('%')
set(gca,'xtick',[])

axes(ax_plot{6})
h=hist(neg_bins,dis_bins);
bar(dis_bins,(h./sum(h))*100,'r')
title(sprintf('Distribution of significant bins\n that are lower than shuffle'))
ylabel('%')
set(gca,'xtick',[])
%% plot z score mat and averages:
mat_dir=cell2mat({per_field_z_score{1}', per_field_z_score{2}'})';

% 1. plot sorted according to pick:
%---------------------------------------
mat_for_sorting=mat_dir(:,1+bins_to_remove:end-bins_to_remove);

%sort:
[~,ind_max]=max(mat_for_sorting,[],2);
[~,ind_sorted]=sort(ind_max);
sorted_mat=mat_dir(ind_sorted,:);
%plot
axes(ax_plot{1})
c_max=z_score_lim;
sorted_mat(find(sorted_mat>c_max))=c_max;
sorted_mat=sorted_mat-min(sorted_mat(:));
plot_2d_mat(dis_bins,1:length(sorted_mat),sorted_mat,c_max)
xlabel('Distance between bats')
ylabel('fields')
title('Per field z score sorted by peak')
%

%2. plot sorted according to min:
%---------------------------------------
%sort:
[~,ind_max]=min(mat_for_sorting,[],2);
[~,ind_sorted]=sort(ind_max);
sorted_mat=mat_dir(ind_sorted,:);

%plot
axes(ax_plot{2})
c_max=z_score_lim;
sorted_mat(find(sorted_mat>c_max))=c_max;
sorted_mat=sorted_mat-min(sorted_mat(:));

plot_2d_mat(dis_bins,1:length(sorted_mat),sorted_mat,c_max)
xlabel('Distance between bats')
ylabel('fields')
title('Per field z score sorted by trough')


%% save
dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end

fig_name=fullfile(dir_name,'population_per_field');
print(gcf, fig_name, '-dpdf',  '-painters');
%% just hists

figure('units','normalized','outerposition',[0 0 1 1])
paper_size = [24 21]; % ~A4 size
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paper_size]);
set(gcf,'PaperOrientation','landscape');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 1 0 0]); % position on screen...
set(gcf, 'Renderer', 'opengl');

ax_plot{4}=axes('position',[0.08 0.35 0.2 0.2]);%hist
ax_plot{6}=axes('position',[0.58 0.35 0.2 0.2]);%hist

%% hist of significant bins
fsize=14;
pos_bins=dis_bins(cell2mat(pos_signif));
neg_bins=dis_bins(cell2mat(neg_signif));

axes(ax_plot{4})
h=hist(pos_bins,dis_bins);
bar(dis_bins,(h./sum(h))*100,'g')
title(sprintf('Distribution of significant bins\n that are higher than shuffle'),'fontsize',fsize,'fontweight','bold')
ylabel('% of significant bins','fontsize',fsize,'fontweight','bold')
xlabel('Inter-bat distance','fontsize',fsize,'fontweight','bold')
set(gca,'xtick',-40:20:40)

axes(ax_plot{6})
h=hist(neg_bins,dis_bins);
bar(dis_bins,(h./sum(h))*100,'r')
title(sprintf('Distribution of significant bins\n that are lower than shuffle'),'fontsize',fsize,'fontweight','bold')
ylabel('% of significant bins','fontsize',fsize,'fontweight','bold')
xlabel('Inter-bat distance','fontsize',fsize,'fontweight','bold')
set(gca,'xtick',-40:20:40)
%% save

dir_name=['D:\Ayelet\2bat_proj\Analysis\new_code\figures\poster_figures\'];
if ~exist(dir_name)
    mkdir(dir_name)
end

fig_name=fullfile(dir_name,'hist_population_per_field');
print(gcf, fig_name, '-dpdf',  '-painters');