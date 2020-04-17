clear
% parameters:
dir_data='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis\';
dir_info=dir(dir_data);
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
load(fullfile(param_folder,'solo_params.mat'));
y_analysis_folder='D:\Ayelet\2bat_proj\Analysis\new_code\figures\y_analysis\';

figure('units','normalized','outerposition',[0 0 1 1],'color',[1 1 1])
%p = parpool(6);
%shuffle params:
n_shuffles=100;
%cell selection:
SI_threshold=1;
min_n_spike=50;
solo_time_spent_minimum_for_1D_bins=0.1;%need to check!
range=[-40 40];
window=10;
step=2;
bin_start=range(1):step:range(2)-window;
bin_end=range(1)+window:step:range(2);
bin_center=[bin_start+bin_end]./2;
bin_dis=range(1):step:range(2);

n_parts=3;

solo_X_n_bins = solo_X_max / 2;
solo_X_bin_size = (solo_X_max-solo_X_min)/solo_X_n_bins;
solo_X_bins_vector=solo_X_min:solo_X_bin_size:solo_X_max;
solo_X_bins_vector_of_centers=solo_X_bins_vector(1:end-1)+solo_X_bin_size/2;

%% 1. load all data:
load(fullfile(param_folder,'solo_params.mat'));
solo_X_n_bins = solo_X_max*2;
solo_X_bin_size = (solo_X_max-solo_X_min)/solo_X_n_bins;
solo_X_bins_vector=solo_X_min:solo_X_bin_size:solo_X_max;
solo_X_bins_vector_of_centers=solo_X_bins_vector(1:end-1)+solo_X_bin_size/2;
%figure prop:
width=0.2;
hight=0.15;
ver_dis=0.4;
hor_dis=0.5;

x_pos(1)=0.05;
x_pos(2)=x_pos(1)+hor_dis;

y_pos(1)=0.6;
y_pos(2)=y_pos(1)-ver_dis;

cell_count=0;

PV_co_dir_bin=[];PV_solo_dir=[];
%%
count=0;

for cell_i=3:length({dir_info.name})-1
    for dir_i=1:2
        
        % for cell_i=3:15
        
        % load data:
        load(fullfile(dir_data,dir_info(cell_i).name))
        bat=cell_co_solo_initial_analysis.exp_data.bat;
        day=cell_co_solo_initial_analysis.exp_data.day;
        cell_num=cell_co_solo_initial_analysis.exp_data.cell_num;
        
        
        num_spike_during_flight=sum(~isnan(cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec(:)))+sum(~isnan(cell_co_solo_initial_analysis.co(dir_i).spikes.ts_usec(:)));
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
        solo_x_pos_bsp=cell_co_solo_initial_analysis.solo(dir_i).bsp.x_pos';solo_x_pos_bsp=solo_x_pos_bsp(:);
        solo_y_pos_bsp=cell_co_solo_initial_analysis.solo(dir_i).bsp.y_pos';solo_y_pos_bsp=solo_y_pos_bsp(:);
        co_x_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.x_pos';co_x_pos_bsp=co_x_pos_bsp(:);
        co_y_pos_bsp=cell_co_solo_initial_analysis.co(dir_i).bsp.y_pos';co_y_pos_bsp=co_y_pos_bsp(:);
     
        
        solo_x_pos_spike=cell_co_solo_initial_analysis.solo(dir_i).spikes.x_pos';solo_x_pos_spike=solo_x_pos_spike(:);
        solo_y_pos_spike=cell_co_solo_initial_analysis.solo(dir_i).spikes.y_pos';solo_y_pos_spike=solo_y_pos_spike(:);
        solo_flight_spike_ts=cell_co_solo_initial_analysis.solo(dir_i).spikes.ts_usec';solo_flight_spike_ts=solo_flight_spike_ts(:);
        %remove nans
        solo_x_pos_bsp=solo_x_pos_bsp(~isnan(solo_x_pos_bsp));
        solo_y_pos_bsp=solo_y_pos_bsp(~isnan(solo_y_pos_bsp));
        co_x_pos_bsp=co_x_pos_bsp(~isnan(co_x_pos_bsp));
        co_y_pos_bsp=co_y_pos_bsp(~isnan(co_y_pos_bsp));
        solo_x_pos_spike=solo_x_pos_spike(~isnan(solo_x_pos_spike));
        solo_y_pos_spike=solo_y_pos_spike(~isnan(solo_y_pos_spike));
        solo_flight_spike_ts=solo_flight_spike_ts(~isnan(solo_flight_spike_ts));
        
        x_pos_per_y_bin_percntile=cell(n_parts,1);
        y_pos_per_y_bin_percntile=cell(n_parts,1);
        x_spike_pos_per_y_bin_percntile=cell(n_parts,1);
        x_pos_per_y_bin_pos=cell(n_parts,1);
        y_pos_per_y_bin_pos=cell(n_parts,1);
        x_spike_pos_per_y_bin_pos=cell(n_parts,1);
        ts_spike_ts_per_y_bin_pos=cell(n_parts,1);
        ts_spike_per_y_bin_percntile=cell(n_parts,1);
        co_10_perc=nan*zeros(size(solo_X_bins_vector_of_centers));
        co_90_perc=nan*zeros(size(solo_X_bins_vector_of_centers));
        solo_10_perc=nan*zeros(size(solo_X_bins_vector_of_centers));
        solo_90_perc=nan*zeros(size(solo_X_bins_vector_of_centers));
        intersect_by_union_cell=nan*zeros(size(solo_X_bins_vector_of_centers));
        
        [~, ~, ~, x_tuning_all_y_data(:,cell_count,dir_i), ~,~] ...
            = fn_compute_generic_1D_tuning_new_smooth ...
            (solo_x_pos_bsp, solo_x_pos_spike, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
        % intersect/union of co and solo y pos
        
        
        % divide data to 3 parts of y per x position
        for bin_i=1:length(solo_X_bins_vector_of_centers)
            relevant_bsp_ind=find(solo_x_pos_bsp>=solo_X_bins_vector(bin_i)& solo_x_pos_bsp<solo_X_bins_vector(bin_i+1));
            if ~isempty(relevant_bsp_ind)
                x_bsp_bin_pos=solo_x_pos_bsp(relevant_bsp_ind);
                relevant_spike_ind=find(solo_x_pos_spike>=solo_X_bins_vector(bin_i)& solo_x_pos_spike<solo_X_bins_vector(bin_i+1));
                x_spike_bin_pos=solo_x_pos_spike(relevant_spike_ind);
                y_spike_bin_pos=solo_y_pos_spike(relevant_spike_ind);
                ts_spike_bin=solo_flight_spike_ts(relevant_spike_ind);
                
                y_bsp_pos_per_x_bin=solo_y_pos_bsp(relevant_bsp_ind);
                solo_10_perc(bin_i)=prctile(y_bsp_pos_per_x_bin,10);
                solo_90_perc(bin_i)=prctile(y_bsp_pos_per_x_bin,90);
                % find relevant y for co:
                co_relevant_bsp_ind=find(co_x_pos_bsp>=solo_X_bins_vector(bin_i)& co_x_pos_bsp<solo_X_bins_vector(bin_i+1));
                co_y_bsp_pos_per_x_bin=co_y_pos_bsp(co_relevant_bsp_ind);
                co_10_perc(bin_i)=prctile(co_y_bsp_pos_per_x_bin,10);
                co_90_perc(bin_i)=prctile(co_y_bsp_pos_per_x_bin,90);
                
                if sum(isnan([co_10_perc(bin_i),solo_10_perc(bin_i),co_90_perc(bin_i),solo_90_perc(bin_i)]))==0
                range_intersect=[max([co_10_perc(bin_i),solo_10_perc(bin_i)]), min([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                range_union=[min([co_10_perc(bin_i),solo_10_perc(bin_i)]), max([co_90_perc(bin_i),solo_90_perc(bin_i)])];
                intersect_by_union_cell(bin_i)=abs(range_intersect(2)-range_intersect(1))/abs(range_union(2)-range_union(1));
                %check if there is no intersection
                if solo_10_perc(bin_i)>co_90_perc(bin_i) | co_10_perc(bin_i)>solo_90_perc(bin_i)
                    intersect_by_union_cell(bin_i)=0;
                end
                else
                    intersect_by_union_cell(bin_i)=nan;
                end
                range_of_y=max(y_bsp_pos_per_x_bin)-min(y_bsp_pos_per_x_bin);
                y_vec_percntile=[min(y_bsp_pos_per_x_bin),prctile(y_bsp_pos_per_x_bin,100/n_parts),prctile(y_bsp_pos_per_x_bin,2*(100/n_parts)),max(y_bsp_pos_per_x_bin)];
                y_vec_pos=[min(y_bsp_pos_per_x_bin),min(y_bsp_pos_per_x_bin)+range_of_y/n_parts,min(y_bsp_pos_per_x_bin)+2*range_of_y/(n_parts),max(y_bsp_pos_per_x_bin)];
                for y_bin=1:n_parts
                    x_pos_per_y_bin_percntile{y_bin}= [x_pos_per_y_bin_percntile{y_bin}; x_bsp_bin_pos(find(y_bsp_pos_per_x_bin>=y_vec_percntile(y_bin) & y_bsp_pos_per_x_bin<y_vec_percntile(y_bin+1)))];
                    y_pos_per_y_bin_percntile{y_bin}= [y_pos_per_y_bin_percntile{y_bin}; y_bsp_pos_per_x_bin(find(y_bsp_pos_per_x_bin>=y_vec_percntile(y_bin) & y_bsp_pos_per_x_bin<y_vec_percntile(y_bin+1)))];

                    x_spike_pos_per_y_bin_percntile{y_bin}= [x_spike_pos_per_y_bin_percntile{y_bin}; x_spike_bin_pos(find(y_spike_bin_pos>=y_vec_percntile(y_bin) & y_spike_bin_pos<y_vec_percntile(y_bin+1)))];
                    ts_spike_per_y_bin_percntile{y_bin}=[ts_spike_per_y_bin_percntile{y_bin};ts_spike_bin(find(y_spike_bin_pos>=y_vec_percntile(y_bin) & y_spike_bin_pos<y_vec_percntile(y_bin+1)))];
                    
                    x_pos_per_y_bin_pos{y_bin}= [x_pos_per_y_bin_pos{y_bin}; x_bsp_bin_pos(find(y_bsp_pos_per_x_bin>=y_vec_pos(y_bin) & y_bsp_pos_per_x_bin<y_vec_pos(y_bin+1)))];
                    y_pos_per_y_bin_pos{y_bin}= [y_pos_per_y_bin_pos{y_bin}; y_bsp_pos_per_x_bin(find(y_bsp_pos_per_x_bin>=y_vec_pos(y_bin) & y_bsp_pos_per_x_bin<y_vec_pos(y_bin+1)))];

                    x_spike_pos_per_y_bin_pos{y_bin}= [x_spike_pos_per_y_bin_pos{y_bin}; x_spike_bin_pos(find(y_spike_bin_pos>=y_vec_pos(y_bin) & y_spike_bin_pos<y_vec_pos(y_bin+1)))];
                    ts_spike_ts_per_y_bin_pos{y_bin}=[ts_spike_ts_per_y_bin_pos{y_bin};ts_spike_bin(find(y_spike_bin_pos>=y_vec_pos(y_bin) & y_spike_bin_pos<y_vec_pos(y_bin+1)))];
                    
                end
            end
        end
        intersect_by_union(cell_count,dir_i)=nanmean(intersect_by_union_cell);
        %compute psth per y bin:
        for y_bin=1:n_parts
            x_bsp=x_pos_per_y_bin_percntile{y_bin};
            x_spike=x_spike_pos_per_y_bin_percntile{y_bin};
            
            [~, ~, ~, percntile_tuning(:,cell_count,y_bin,dir_i), ~,~] ...
                = fn_compute_generic_1D_tuning_new_smooth ...
                (x_bsp, x_spike, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
            
            x_bsp=x_pos_per_y_bin_pos{y_bin};
            x_spike=x_spike_pos_per_y_bin_pos{y_bin};
            
            [~, ~, ~, pos_tuning(:,cell_count,y_bin,dir_i), ~,~] ...
                = fn_compute_generic_1D_tuning_new_smooth ...
                (x_bsp, x_spike, solo_X_bins_vector_of_centers, solo_time_spent_minimum_for_1D_bins, frames_per_second, 0,0,0);
            
            
        end
        %compute corrs:
        %precntile:
        corr_mat=corr(squeeze(percntile_tuning(:,cell_count,:,dir_i)), 'Rows', 'pairwise');
        corrs_percntile_bins=tril(corr_mat,-1);
        corrs_percntile_bins=corrs_percntile_bins(:);
        corrs_percntile_bins(find(corrs_percntile_bins==0))=[];
        
        corr_percntile_all=[corr_mat, corr(squeeze(percntile_tuning(:,cell_count,:,dir_i)),squeeze(x_tuning_all_y_data(:,cell_count,dir_i)), 'Rows', 'pairwise')];
        
        %pos range:
        corr_mat=corr(squeeze(pos_tuning(:,cell_count,:,dir_i)), 'Rows', 'pairwise');
        corrs_percntile_bins=tril(corr_mat,-1);
        corrs_percntile_bins=corrs_percntile_bins(:);
        corrs_percntile_bins(find(corrs_percntile_bins==0))=[];
        
        corr_pos_all=[corr_mat, corr(squeeze(pos_tuning(:,cell_count,:,dir_i)),x_tuning_all_y_data(:,cell_count,dir_i), 'Rows', 'pairwise')];
        %% plots
        %title
        if dir_i==1
            % write cell's ID
            str = sprintf('Cell #%d, Date: %d, Bat %d',cell_num,day,bat) ;
            annotation('textbox',[.06 .88 .3 .1],'string',str,'EdgeColor','none','fontsize',15)
        end
        % percentile:
        axes('position',[x_pos(dir_i) y_pos(1),width,hight])
        plot(solo_X_bins_vector_of_centers,x_tuning_all_y_data(:,cell_count,dir_i),'k')
        hold on
        plot(solo_X_bins_vector_of_centers,squeeze(percntile_tuning(:,cell_count,:,dir_i)))
        ylabel('Hz')
        title('divide by perctile of y pos')
        xlim([0 135])
        
        axes('position',[x_pos(dir_i) y_pos(1)-hight,width,hight])
        for y_bin=1:n_parts
            ts=(ts_spike_per_y_bin_percntile{y_bin}-min(solo_flight_spike_ts))/1e6/60;
            plot(x_spike_pos_per_y_bin_percntile{y_bin},ts_spike_per_y_bin_percntile{y_bin},'.'); hold on
        end
        ylabel('Time (min)')
        xlim([0 135])
        
        axes('position',[x_pos(dir_i)+width*1.3 y_pos(1)-hight,width*0.8,hight])
        x = repmat(1:size(corr_percntile_all,2),size(corr_percntile_all,1),1); % generate x-coordinates
        y = repmat((1:size(corr_percntile_all,1))',1,size(corr_percntile_all,2)); % generate y-coordinates
        % Generate Labels
        t = num2cell(corr_percntile_all); % extact values into cells
        t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
        % Draw Image and Label Pixels
        imagesc(corr_percntile_all)
        text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
        set(gca,'xtick',1:size(corr_percntile_all,2),'xticklabel',{'part 1','part 2','part 3','all'})
        set(gca,'ytick',1:size(corr_percntile_all,2),'yticklabel',{'part 1','part 2','part 3'}   )
        colorbar
        title('corr values between parts percentile')
        
        axes('position',[x_pos(dir_i)+width*1.3 y_pos(1)+0.05,width*0.8,hight*0.8])
        plot(x_pos_per_y_bin_percntile{1},y_pos_per_y_bin_percntile{1},'.');hold on;
        plot(x_pos_per_y_bin_percntile{2},y_pos_per_y_bin_percntile{2},'.');hold on;
        plot(x_pos_per_y_bin_percntile{3},y_pos_per_y_bin_percntile{3},'.');hold on;
        xlim([0 135])
        xlabel('x pos (m)')
        ylabel('y pos (m)')
        
        %divide by pos
        axes('position',[x_pos(dir_i), y_pos(2),width,hight])
        plot(solo_X_bins_vector_of_centers,x_tuning_all_y_data(:,cell_count,dir_i),'k')
        hold on
        plot(solo_X_bins_vector_of_centers,squeeze(pos_tuning(:,cell_count,:,dir_i)))
        ylabel('Hz')
        title('divide by third of y pos range')
        xlim([0 135])
        
        axes('position',[x_pos(dir_i),y_pos(2)-hight,width,hight])
        for y_bin=1:n_parts
            ts=(ts_spike_ts_per_y_bin_pos{y_bin}-min(solo_flight_spike_ts))/1e6/60;
            plot(x_spike_pos_per_y_bin_pos{y_bin},ts_spike_ts_per_y_bin_pos{y_bin},'.'); hold on
        end
        ylabel('Time (min)')
        xlim([0 135])
        
        
        axes('position',[x_pos(dir_i)+width*1.3 y_pos(2)-hight,width*0.8,hight])
        x = repmat(1:size(corr_pos_all,2),size(corr_pos_all,1),1); % generate x-coordinates
        y = repmat((1:size(corr_pos_all,1))',1,size(corr_pos_all,2)); % generate y-coordinates
        % Generate Labels
        t = num2cell(corr_pos_all); % extact values into cells
        t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
        % Draw Image and Label Pixels
        imagesc(corr_pos_all)
        text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
        set(gca,'xtick',1:size(corr_pos_all,2),'xticklabel',{'part 1','part 2','part 3','all'})
        set(gca,'ytick',1:size(corr_pos_all,2),'yticklabel',{'part 1','part 2','part 3'}   )
        colorbar
        title('corr values between parts pos range')
        
        axes('position',[x_pos(dir_i)+width*1.3 y_pos(2)+0.05,width*0.8,hight*0.8])
        plot(x_pos_per_y_bin_pos{1},y_pos_per_y_bin_pos{1},'.');hold on;
        plot(x_pos_per_y_bin_pos{2},y_pos_per_y_bin_pos{2},'.');hold on;
        plot(x_pos_per_y_bin_pos{3},y_pos_per_y_bin_pos{3},'.');hold on;
        xlim([0 135])
        xlabel('x pos (m)')
        ylabel('y pos (m)')
        
        axes('position',[x_pos(dir_i)+width*1.3 0.8,width*0.8,hight*0.8])
        plot(solo_X_bins_vector_of_centers,[co_10_perc;co_90_perc],'m');hold on;
       plot(solo_X_bins_vector_of_centers,[solo_10_perc;solo_90_perc],'k')
                      ylabel(['intersect/union',num2str(intersect_by_union(cell_count,dir_i))])

      axes('position',[x_pos(dir_i)+width*1.3 0.8+hight*0.8,width*0.8,hight*0.5])
        plot(solo_X_bins_vector_of_centers,intersect_by_union_cell)

    end
         %save figure
        
        if ~exist(y_analysis_folder)
            mkdir(y_analysis_folder)
        end
        fig_name=fullfile(y_analysis_folder,['cell_',num2str(cell_num),'_day_',num2str(day),'_bat_',num2str(bat),'.png']);
        saveas(gcf,fig_name)
        clf
end


%% population
for dir_i=1:2
    [percntile_cell_corr_1_2(:,dir_i) percntile_cell_corr_mean_1_2]=corr_mat_cells_by_pos(percntile_tuning(:,:,1,dir_i),percntile_tuning(:,:,2,dir_i),1);
    [percntile_cell_corr_1_3(:,dir_i) percntile_cell_corr_mean_1_3]=corr_mat_cells_by_pos(percntile_tuning(:,:,1),percntile_tuning(:,:,3,dir_i),1);
    [percntile_cell_corr_3_2(:,dir_i) percntile_cell_corr_mean_3_2]=corr_mat_cells_by_pos(percntile_tuning(:,:,3),percntile_tuning(:,:,2,dir_i),1);
    
    
    [percntile_cell_corr_1_all(:,dir_i) percntile_cell_corr_mean_1_all]=corr_mat_cells_by_pos(percntile_tuning(:,:,1,dir_i),x_tuning_all_y_data(:,:,dir_i),1);
    [percntile_cell_corr_2_all(:,dir_i) percntile_cell_corr_mean_1_all]=corr_mat_cells_by_pos(percntile_tuning(:,:,1,dir_i),x_tuning_all_y_data(:,:,dir_i),1);
    [percntile_cell_corr_3_all(:,dir_i) percntile_cell_corr_mean_3_all]=corr_mat_cells_by_pos(percntile_tuning(:,:,3,dir_i),x_tuning_all_y_data(:,:,dir_i),1);
    
    
    [pos_cell_corr_1_2(:,dir_i) pos_cell_corr_mean_1_2]=corr_mat_cells_by_pos(pos_tuning(:,:,1,dir_i),pos_tuning(:,:,2,dir_i),1);
    [pos_cell_corr_1_3(:,dir_i) pos_cell_corr_mean_1_3]=corr_mat_cells_by_pos(pos_tuning(:,:,1,dir_i),pos_tuning(:,:,3,dir_i),1);
    [pos_cell_corr_3_2(:,dir_i) pos_cell_corr_mean_3_2]=corr_mat_cells_by_pos(pos_tuning(:,:,3,dir_i),pos_tuning(:,:,2,dir_i),1);
    

    [pos_cell_corr_1_all(:,dir_i) pos_cell_corr_mean_1_all]=corr_mat_cells_by_pos(pos_tuning(:,:,1,dir_i),x_tuning_all_y_data(:,:,dir_i),1);
    [pos_cell_corr_2_all(:,dir_i) pos_cell_corr_mean_1_all]=corr_mat_cells_by_pos(pos_tuning(:,:,1,dir_i),x_tuning_all_y_data(:,:,dir_i),1);
    [pos_cell_corr_3_all(:,dir_i) pos_cell_corr_mean_3_all]=corr_mat_cells_by_pos(pos_tuning(:,:,3,dir_i),x_tuning_all_y_data(:,:,dir_i),1);
    
end

%% plot
%figure prop:
width=0.15;
hight=0.1;
ver_dis=0.15;
hor_dis=0.2;


x_pos(1)=0.05;
x_pos(2)=x_pos(1)+hor_dis;
x_pos(3)=x_pos(2)+hor_dis;
x_pos(4)=x_pos(3)+hor_dis;

y_pos(1)=0.7;
y_pos(2)=y_pos(1)-ver_dis;
y_pos(3)=y_pos(2)-ver_dis;
y_pos(4)=y_pos(3)-ver_dis;


%precentile
axes('position',[x_pos(1) y_pos(1),width,hight])
[h x]=hist(percntile_cell_corr_1_2(:))
bar(x,h/sum(h))
title(['corr part 1+2 percentile mean=',num2str(nanmean(percntile_cell_corr_1_2(:)))])
xlim([-1 1])
axes('position',[x_pos(2) y_pos(1),width,hight])
[h x]=hist(percntile_cell_corr_1_3(:))
bar(x,h/sum(h))
title(['corr part 1+3 percentile mean=',num2str(nanmean(percntile_cell_corr_1_3(:)))])
xlim([-1 1])

axes('position',[x_pos(3) y_pos(1),width,hight])
[h x]=hist(percntile_cell_corr_3_2(:))

bar(x,h/sum(h))
title(['corr part 2+3 percentile mean=',num2str(nanmean(percntile_cell_corr_3_2(:)))])
xlim([-1 1])

axes('position',[x_pos(1) y_pos(2),width,hight])
[h x]=hist(percntile_cell_corr_1_all(:))
bar(x,h/sum(h))
title(['corr part 1+all percentile mean=',num2str(nanmean(percntile_cell_corr_1_all(:)))])
xlim([-1 1])

axes('position',[x_pos(2) y_pos(2),width,hight])
[h x]=hist(percntile_cell_corr_2_all(:))
bar(x,h/sum(h))
title(['corr part 2+all percentile mean=',num2str(nanmean(percntile_cell_corr_2_all(:)))])
xlim([-1 1])

axes('position',[x_pos(3) y_pos(2),width,hight])
[h x]=hist(percntile_cell_corr_3_all(:))
bar(x,h/sum(h))
title(['corr part 3+all percentile mean=',num2str(nanmean(percntile_cell_corr_3_all(:)))])
xlim([-1 1])

%pos range
axes('position',[x_pos(1) y_pos(3),width,hight])
[h x]=hist(pos_cell_corr_1_2(:))
bar(x,h/sum(h))
title(['corr part 1+2 pos range mean=',num2str(nanmean(pos_cell_corr_1_2(:)))])
xlim([-1 1])

axes('position',[x_pos(2) y_pos(3),width,hight])
[h x]=hist(pos_cell_corr_1_3(:))
bar(x,h/sum(h))
title(['corr part 1+3 pos range mean=',num2str(nanmean(pos_cell_corr_1_3(:)))])
xlim([-1 1])

axes('position',[x_pos(3) y_pos(3),width,hight])
[h x]=hist(pos_cell_corr_3_2(:))
bar(x,h/sum(h))
title(['corr part 2+3 pos range mean=',num2str(nanmean(pos_cell_corr_3_2(:)))])
xlim([-1 1])

axes('position',[x_pos(1) y_pos(4),width,hight])
[h x]=hist(pos_cell_corr_1_all(:))
bar(x,h/sum(h))
title(['corr part 1+all pos range mean=',num2str(nanmean(pos_cell_corr_1_all(:)))])
xlim([-1 1])

axes('position',[x_pos(2) y_pos(4),width,hight])
[h x]=hist(pos_cell_corr_2_all(:))
bar(x,h/sum(h))
title(['corr part 2+all pos range mean=',num2str(nanmean(pos_cell_corr_2_all(:)))])
xlim([-1 1])

axes('position',[x_pos(3) y_pos(4),width,hight])
[h x]=hist(pos_cell_corr_3_all(:))
bar(x,h/sum(h))
title(['corr part 3+all pos range mean=',num2str(nanmean(pos_cell_corr_3_all(:)))])
xlim([-1 1])



%intersect/union:
intersect_by_union(intersect_by_union==0)=nan;
axes('position',[x_pos(3)+width*1.1 y_pos(4),width,hight])
[h x]=hist(intersect_by_union)
bar(x,h/sum(h))
title(['intersect/union',num2str(nanmean(intersect_by_union(:)))])
xlim([0 1])






fig_name=fullfile(y_analysis_folder,['population_hists.png']);
        saveas(gcf,fig_name)
        clf
