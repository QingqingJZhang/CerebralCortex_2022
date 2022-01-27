function [on_matrix, off_matrix]=get_BOLD_traces(BOLD,light_onset,light_offset,on_id,off_id)
% for the fixed stimulation, show the response from -5 to 24 s; 
if isempty(on_id)
    on_matrix = [];
else
    fMRI_on_matrix=mean(BOLD(on_id,:),1,'omitnan');    
    on_matrix = zeros(length(light_onset), 30);    
    for i = 1:length(light_onset)
        on_matrix(i,:) = fMRI_on_matrix(light_onset(i)-5:light_onset(i)+24);
    end
end

if isempty(off_id)
    off_matrix = [];
else
    fMRI_off_matrix = mean(BOLD(off_id,:),1,'omitnan');    
    off_matrix = zeros(length(light_offset), 30);    
    for i = 1:length(light_offset)
        off_matrix(i,:) = fMRI_off_matrix(light_offset(i)-5:light_offset(i)+24);
    end
end

end





