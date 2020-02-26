function plot_2d_mat(x_bins,y_bins,mat,c_max)
%mat=mat-min(mat(:));
hhh=imagesc(x_bins,y_bins,mat);
colormap_map=jet;
colormap_map(64,:) = [1 1 1] ; % Set the highest-value pixel to White = [1 1 1] : I will use White as the Color of NaN pixels
%if isempty(c_max)
cdata_mat = mat / max(mat(:)) * ...
    ( size(colormap_map,1) - 1 ) ; % Scale the colors in 'cdata' to ( 1 : length(colormap) - 1 ), so that the highest value will be reserved to NaN's
% else
%     
% %set all that above c_max to c_max
% mat(find(mat>c_max))=size(colormap_map,1) - 1;
%  cdata_mat = mat / max(mat(:)) * ...
%     ( size(colormap_map,1) - 1 ) ; % Scale the colors in 'cdata' to ( 1 : length(colormap) - 1 ), so that the highest value will be reserved to NaN's
% end
idx_isNaN_Field_X_Y = find( isnan( mat  ) ); % Find the indexes of NaN bins
cdata_mat( idx_isNaN_Field_X_Y ) = size(colormap_map,1) ; % Replace NaN values with the highest values (they will be colored white)
colormap( colormap_map );
caxis([0 65]); % Scale the lowest value (deep blue) to 0
set(hhh, 'cdatamapping', 'direct', 'cdata', cdata_mat );
axis xy