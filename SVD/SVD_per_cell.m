% SVD analysis
sigma_a=1.5;
hsize=5*round(sigma_a)+1;
%------------------------------------------------------------------------
% run over cells:

%% structs' parameters
cell_co_solo_initial_analysis_struct_folder='D:\Ayelet\2bat_proj\Analysis\new_code\analysis_structs\co_solo_initial_analysis';
files = dir(cell_co_solo_initial_analysis_struct_folder);
behavior_struct_names = {files.name};
is_dir=[files.isdir];
behavior_struct_names=behavior_struct_names(~is_dir);

figure;
co_fig_folder_name='D:\Ayelet\2bat_proj\Analysis\new_code\figures\SVD\';

%% plot for each cell
for ii_cell =1:length(behavior_struct_names)
    ii_cell
    %try
    signif=0;
    %% load data
    
    struct_name =behavior_struct_names{ii_cell};
    file_name = fullfile(cell_co_solo_initial_analysis_struct_folder,struct_name);
    load(file_name);
    behavior_struct=cell_co_solo_initial_analysis;
    bat=behavior_struct.exp_data.bat;
    day=behavior_struct.exp_data.day;
    if isnumeric(day)
        day=num2str(day);
    end
    cell_num=behavior_struct.exp_data.cell_num;
    % read cell mat
    for dir_i=1:2
        x=cell_co_solo_initial_analysis.co(dir_i).firing_rate.field_density_smoothed_XY_with_NaN;
        x(sum(isnan(x),2)==size(x,2),:)=[];
        if isempty(x)
            continue
        end
        % guess initial filling of nans:
        nonnan_ind = ~isnan(x);
        nan_ind=isnan(x);
        vvalues = x(nonnan_ind);
        %bsxfun compute the distance from every nan to all non-nan.
        %min then find which of the distance is the smallest
        %         [~, mindistcol] = min(abs(bsxfun(@minus, (1:numel(x))', find(nonnan_ind)')));
                 x_no_nan=x;
        %         x_no_nan(nan_ind)=x(mindistcol(find(nan_ind)));
        
        STATS = regionprops(nan_ind,'Centroid','PixelIdxList','ConvexHull');
        for ii_field=1:1:length(STATS)
%             PixelIdxList(ii_field) =STATS(ii_field).PixelIdxList;
%             XY_field_centroid(ii_field,:) = STATS(ii_field).Centroid;
%             XY_field(ii_field).ConvexHull = STATS(ii_field).ConvexHull;
            p = polyshape( STATS(ii_field).ConvexHull); %create a polygon from x y coordinates
            [cx, cy] = centroid(p); %find the centroid
            factor = 0.5;
            %m = scale(p, factor, [cx cy]) %scale p by factor 0.5 w.r.t a point [cx cy]. m's area is 25% of p's area.
            q = polybuffer(p,+1);
            %x_limits=round([min(q.Vertices(:,1)) max(q.Vertices(:,1))]);
            %y_limits=round([min(q.Vertices(:,2)) max(q.Vertices(:,2))]);
            %x_limits(x_limits<=0)=1;x_limits(x_limits>size(x,2))=size(x,2);
            %y_limits(y_limits<=0)=1;y_limits(y_limits>size(x,1))=size(x,1);
            %
            [xx,yy]=meshgrid(1:size(x,2),1:size(x,1));
            in = inpolygon(xx(:),yy(:),q.Vertices(:,1),q.Vertices(:,2));
            ind=sub2ind(size(x),yy(in),xx(in));%FIX!!!!!!!!
            x_no_nan(STATS(ii_field).PixelIdxList)=nanmean(x_no_nan(ind));
        end
        
        % comute SVD
        try
        [U,S,V] = svd(x_no_nan);
        s=diag(S);
        alpha_svd=1-s(1).^2./sum(s.^2);
        % shuffle
        ego_tuning=nanmean(x_no_nan);
        allo_tuning=nanmean(x_no_nan,2);
        ego_by_allo_map=ego_tuning.*allo_tuning;
        ego_by_allo_map=ego_by_allo_map./max(ego_by_allo_map(:));
        mean_fr=mean(x_no_nan(:));
        ego_by_allo_map=ego_by_allo_map*mean_fr;
        for shuffle_i=1:1000
            r=poissrnd(ego_by_allo_map);
            gaussian_kernel = fspecial('gaussian',hsize,sigma_a);
            
            % Smoothing = convolve with gaussian kernel:
            r_smoothed = imfilter(r,gaussian_kernel);
            
            [U,S,V] = svd(r_smoothed);
            s=diag(S);
            alpha_svd_shuffle(shuffle_i)=1-s(1).^2./sum(s.^2);
        end
        %%
        subplot(2,2,1)
        imagesc(x_no_nan)
        title('real map')
        subplot(2,2,2)
        imagesc(ego_by_allo_map)
        title('separable map')
        subplot(2,2,3)
        hist(alpha_svd_shuffle); hold on; plot([alpha_svd alpha_svd],[0 300],'r')
        alpha_svd_all=[alpha_svd,alpha_svd_shuffle];
        alpha_svd_all_sort=sort(alpha_svd_all);
        p=1-find(alpha_svd_all_sort==alpha_svd)./length(alpha_svd_all_sort);
        title(sprintf('alpha SVD=%.3f , p=%.3f',alpha_svd,p))
        %%
        fig_name=fullfile(co_fig_folder_name,['cell_',num2str(cell_num),'_dir_',num2str(dir_i),'_day_',day,'_bat_',num2str(bat),'.png']);
        saveas(gcf,fig_name)
        clf
        catch
        end
    end
end
% while for completing the
% x=cell_co_solo_initial_analysis.co(1).firing_rate.field_density_smoothed_XY_with_NaN;
% x(sum(isnan(x),2)==26,:)=[];
% find_nan_mat=zeros(size(x));
% find_nan_mat(isnan(x))=1;
% % % Find individual fields
% STATS = regionprops(find_nan_mat,'Centroid','PixelIdxList','ConvexHull');
% for ii_field=1:1:length(STATS)
%     XY_field_size(ii_field) = length(STATS(ii_field).PixelIdxList);
%     XY_field_centroid(ii_field,:) = STATS(ii_field).Centroid;
%     XY_field(ii_field).ConvexHull = STATS(ii_field).ConvexHull;
%
% end;
% % inpolygon

%