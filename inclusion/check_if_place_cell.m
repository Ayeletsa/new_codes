function inclusion=check_if_place_cell(inclusion,solo,ii_dir,shuffling_struct,param_file_name)
load(param_file_name)
if inclusion(ii_dir).valid_cell==0
    inclusion(ii_dir).place_cell=0;
else

if isempty(solo(ii_dir).fields)
    has_valid_field = 0;
else
    has_valid_field = 1;
end
% SI > thr
SI_thr_signif = solo(ii_dir).SI   > SI_thr_solo;
% SI > shuffle
SI_shuffle_signif = solo(ii_dir).SI > ...
    prctile([shuffling_struct(ii_dir).FE_PSTH_shuffle.FE_PSTH.SI_bits_spike],SI_thr_shuffle_solo);
% valid cell
valid_cell=inclusion(ii_dir).valid_cell;
%pyramidal:
pyr=inclusion(ii_dir).pyr;
%% apply all conditions
TF = true;
TF = TF & SI_thr_signif;
TF = TF & SI_shuffle_signif;
TF = TF & has_valid_field;
TF = TF & pyr; 
TF = TF & valid_cell; 

inclusion(ii_dir).place_cell=TF;
end