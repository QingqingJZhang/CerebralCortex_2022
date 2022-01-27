
function [on_matrix, off_matrix]=get_BOLD_traces_ran(BOLD,light_onset,light_offset,on_id,off_id)

if isempty(on_id)
    on_matrix = [];
else
    fMRI_on_matrix=mean(BOLD(on_id,:),1,'omitnan');    
    on_matrix = {};    
    for i = 1:length(light_onset)
        if i >length(light_offset)
            on_matrix{i} = fMRI_on_matrix(light_onset(i)-5:end);
        else
            on_matrix{i} = fMRI_on_matrix(light_onset(i)-5:light_offset(i));
        end
    end
end

if isempty(off_id)
    off_matrix = [];
else
    fMRI_off_matrix = mean(BOLD(off_id,:),1,'omitnan');    
    off_matrix = {};    
    for i = 1:length(light_offset)
        if i+1>length(light_onset)
            off_matrix{i} = fMRI_off_matrix(light_offset(i)-5:end);
        else
            off_matrix{i} = fMRI_off_matrix(light_offset(i)-5:light_onset(i+1));
        end
    end
end

end





