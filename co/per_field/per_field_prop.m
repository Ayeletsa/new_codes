function per_field=per_field_prop(per_field,field_i,r,p)


per_field(field_i).sparsity=(nanmean(r).^2)/nanmean(r.^2);
std_r=nanstd(r);
mean_r=nanmean(r);
per_field(field_i).cv=std_r/mean_r;
per_field(field_i).modulation_depth=(max(r)-min(r))./(max(r));
per_field(field_i).SI = fn_compute_spatial_info (p,r);
