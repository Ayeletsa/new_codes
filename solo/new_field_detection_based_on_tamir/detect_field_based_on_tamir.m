function fields=detect_field_based_on_tamir(FR_map,FE,prm) 
    
    %% detect peaks
    fields = fields_detect_peaks(FR_map,prm);
    
    %% Field width/edges
    [widths, edges] = fields_calc_width_edges(FR_map, fields, prm.fields.width_href);
    [fields(:).width_href]   = disperse(widths);
    [fields(:).edges_href]   = disperse(edges);

    %% remove overlapping, lower fields
    fields = fields_remove_overlaps(FR_map, fields, prm);
    
    %% limit field size if neighboring fields overlap (in field width and not in ref height for overlap)
    fields = fields_limit_width_by_neighbors(FR_map, fields, prm);
    
    %% get spikes in field
    fields = fields_add_spikes_data(FE, FR_map, fields, prm);
    
    %% remove unstable fields
    fields = fields_remove_unstable(fields, prm);
    
    %% Remove non significant fields based on local shuffling analysis
    fields = fields_find_signif(FE, fields, prm);

    %% sort fields by position
    [~,IX] = sort([fields.loc],'ascend');
    fields = fields(IX);
    
end