function [on_map, off_map, on_off_map] = get_BOLD_correlation_maps(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,nan_index,hrf)

on_matrix = zeros(1,size(BOLD,2));
off_matrix = zeros(1,size(BOLD,2));
on_off_matrix = zeros(1,size(BOLD,2));

on_matrix(adjusted_epi_frames_on) = 1;
off_matrix(adjusted_epi_frames_off) = 1;
on_off_matrix(adjusted_epi_frames_on) = 1; on_off_matrix(adjusted_epi_frames_off) = 1;
off_design = conv(off_matrix,hrf,'full'); off_design = off_design(1:end-length(hrf)+1);
off_design(nan_index) = [];
BOLD(:,nan_index) = [];
off_map = corr(BOLD', off_design');

on_design = conv(on_matrix,hrf,'full'); on_design = on_design(1:end-length(hrf)+1);
on_design(nan_index) = [];
on_map = corr(BOLD', on_design');

on_off_design = conv(on_off_matrix,hrf,'full'); on_off_design = on_off_design(1:end-length(hrf)+1);
on_off_design(nan_index) = [];
on_off_map = corr(BOLD',on_off_design');

end
