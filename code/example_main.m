% 20210922 Qingqing Zhang 

% Install SPM12 (http://www.fil.ion.ucl.ac.uk/spm/). 
% Change "p" to [3 6 1 2 2 0 32] inside the "spm_hrf" function. 
%% 1.1 Loop through the 4-channel data and extract traces from the SC, VA, Cg
clear;


load('YOURPATH/masks/4chan_rSCrVACg_masks.mat');
% load('YOURPATH/4chan_non_visual_mask.mat');
load('YOURPATH/Big_ACC_anatomy.mat');

path = 'YOURPATH/4-channel';
animals = dir(path);


on_SCs_f = []; off_SCs_f=[];
on_SCs_r = {}; off_SCs_r={};
on_Vs_f = []; off_Vs_f=[];
on_Vs_r = {}; off_Vs_r={};
on_Cgs_f = []; off_Cgs_f=[];
on_Cgs_r = {}; off_Cgs_r={};

for i = 1:size(animals,1)    
    id = strfind(animals(i).name, 'V');
    if ~isempty(id)
        animal_id = animals(i).name(id+1:end); 
        cd(fullfile(path,['WV',num2str(animal_id)]));
        animal_folder = cd;        
        scan_dates = dir(cd);
        for d = 1:length(scan_dates)-2
            epi_folder = fullfile(animal_folder, scan_dates(d+2).name); 
            cd(epi_folder);
            EPIs = dir(cd);
            for e = 1:length(EPIs)-2
                cd(fullfile(epi_folder, EPIs(e+2).name));%
                load('trigger.mat'); load('BOLD');
                % Normalize the BOLD signal to resting-state baseline
                baseline_frames = 1:baseline_TTLs;
                voxel_means = mean(BOLD(:,baseline_frames),2,'omitnan');
                voxel_stds = std(BOLD(:,baseline_frames),[],2,'omitnan');
                target_epis = BOLD(:,length(baseline_frames)-9:end);
                target_epis = target_epis - voxel_means;
                target_epis = target_epis./voxel_stds;               
%                 % regress out mean response of 'nonresponding voxels':
%                 global_signal = mean(target_epis(non_visual_mask,:),1,'omitnan');
%                 global_signal = [ones(length(global_signal),1), global_signal'];
%                 for v=1:size(target_epis,1)
%                     [~,~,target_epis(v,:)] = regress(target_epis(v,:)',global_signal);
%                 end

                BOLD = target_epis;                
                on_V_id = on_V_index{i}; off_V_id = off_V_index{i};
                on_Cg_id = on_ACC_index{i}; off_Cg_id = ACC_id;
                on_SC_id = on_SC_index{i}; off_SC_id = off_SC_index{i};
                
                % determine the stimulation paradigm is fixed or randomized:
                if length(unique(diff(adjusted_epi_frames_on)))<=3 % fixed stimualtion
                    [on_SC, off_SC]=get_BOLD_traces(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,on_SC_id,off_SC_id);
                    on_SCs_f = [on_SCs_f;on_SC];
                    off_SCs_f = [off_SCs_f;off_SC];                    
                    
                    [on_V, off_V]=get_BOLD_traces(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,on_V_id,off_V_id);
                    on_Vs_f = [on_Vs_f;on_V];
                    off_Vs_f = [off_Vs_f;off_V];
                    %
                    [on_Cg, off_Cg]=get_BOLD_traces(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,on_Cg_id,off_Cg_id);
                    on_Cgs_f = [on_Cgs_f;on_Cg];
                    off_Cgs_f = [off_Cgs_f;off_Cg];                    
                    
                else %randomized on off
                    [on_SC, off_SC]=get_BOLD_traces_ran(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,on_SC_id,off_SC_id);
                    on_SCs_r = {on_SCs_r{:},on_SC};
                    off_SCs_r = {off_SCs_r{:},off_SC};
                    
                    [on_V, off_V]=get_BOLD_traces_ran(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,on_V_id,off_V_id);
                    on_Vs_r = {on_Vs_r{:},on_V};
                    off_Vs_r = {off_Vs_r{:},off_V};
                
                    [on_Cg, off_Cg]=get_BOLD_traces_ran(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,on_Cg_id,off_Cg_id);
                    on_Cgs_r = {on_Cgs_r{:},on_Cg};
                    off_Cgs_r = {off_Cgs_r{:},off_Cg};  
                
                end                
                clear BOLD adjusted_epi_frames_off adjusted_epi_frames_on 
            end            
        end
    end    
end

%% 1.2 Plot the traces (combine fixed and randomized stimulation):
% show the time course of SC for example:
temp_on = on_SCs_r; temp_off = off_SCs_r; 
len = 25; % plot -5 to 19 s

on_matrix = [];
for i=1:size(temp_on,2)
    for j = 1:size(temp_on{i},2)
        if length(temp_on{i}{j})>len
            on_matrix = [on_matrix; temp_on{i}{j}(1:len)];
        end
    end
end

off_matrix = [];
for i=1:size(temp_off,2)
    for j = 1:size(temp_off{i},2)
        if length(temp_off{i}{j})>len
            off_matrix = [off_matrix; temp_off{i}{j}(1:len)];
        end
    end
end

on_matrix = [on_matrix; on_SCs_f(:,1:25)]; off_matrix = [off_matrix; off_SCs_f(:,1:25)]; 

figure; plot(mean(on_matrix, 'omitnan')); 
figure; plot(mean(off_matrix, 'omitnan')); 

%% 1.3 Loop through the 4-channel data and calcuate the correlation maps: 
clear;
path = 'YOURPATH/4-channel';
animals = dir(path);
hrf = spm_hrf(1);
cepis = 0;
for i = 1:size(animals,1)    
    id = strfind(animals(i).name, 'V');
    if ~isempty(id)
        animal_id = animals(i).name(id+1:end); 
        cd(fullfile(path,['WV',num2str(animal_id)]));
        animal_folder = cd;        
        scan_dates = dir(cd);
        for d = 1:length(scan_dates)-2
            epi_folder = fullfile(animal_folder, scan_dates(d+2).name); 
            cd(epi_folder);
            EPIs = dir(cd);
            for e = 1:length(EPIs)-2
                cd(fullfile(epi_folder, EPIs(e+2).name));%
                load('trigger.mat'); 
                if length(unique(diff(adjusted_epi_frames_on)))<=3 % for fixed visual stimulation: 
                    cepis = cepis + 1;
                load('BOLD');
                % Normalize the BOLD signal to resting-state baseline
                baseline_frames = 1:baseline_TTLs;
                voxel_means = mean(BOLD(:,baseline_frames),2,'omitnan');
                voxel_stds = std(BOLD(:,baseline_frames),[],2,'omitnan');
                target_epis = BOLD(:,length(baseline_frames)-9:end);
                target_epis = target_epis - voxel_means;
                target_epis = target_epis./voxel_stds;%
                BOLD = target_epis;
                % NaN values indicate animal moved at that time point
                nan_index = isnan(BOLD(1,:));                
                [on_map, off_map, on_off_map]=get_BOLD_correlation_maps(BOLD,adjusted_epi_frames_on,adjusted_epi_frames_off,nan_index,hrf);
                on_maps(:,cepis) = atanh(on_map);
                off_maps(:,cepis) = atanh(off_map);
                on_off_maps(:,cepis) = atanh(on_off_map);
                end               
                
                clear BOLD adjusted_epi_frames_off adjusted_epi_frames_on 
            end            
        end
    end    
end

aver_off_maps = mean(off_maps,2,'omitnan');
show_brain_map_highres(aver_off_maps,[],0,0.07,1,4:18,1);cbr=colorbar; set(cbr, 'YTick', [-0.07,0,0.07]); title('Off correlation map');

aver_on_maps = mean(on_maps,2,'omitnan');
show_brain_map_highres(aver_on_maps,[],0,0.07,1,4:18,1);cbr=colorbar; set(cbr, 'YTick', [-0.07,0,0.07]); title('On correlation map');

aver_on_off_maps = mean(on_off_maps,2,'omitnan');
show_brain_map_highres(aver_on_off_maps,[],0,0.07,1,4:18,1);cbr=colorbar; set(cbr, 'YTick', [-0.07,0,0.07]); title('On off correlation map');

%% 2. Linear model of using left-hemisphere ROIs'response fit rSC's response with the 4-channel data:
% "linearM_4chan_fixedStim.mat" contains trial-by-trial ON OFF BOLD responses from the SC (on_SCs_f; off_SCs_f) and 36 left-hemisphere ROIs; 

clear;
load('linearM_4chan_fixedStim'); % use data from fixed stimulation
on_matrix =  on_SCs_f;
off_matrix = off_SCs_f; 
test_on = on_lROIs_f; % ON response from each left-hemisphere ROI
test_off = off_lROIs_f; % OFF response from each left-hemisphere ROI
% remove trials with motions
on_matrix(isnan(on_matrix))=0; on_zero = ~all(on_matrix,2);on_matrix(on_zero,:) = [];
off_matrix(isnan(off_matrix))=0; off_zero = ~all(off_matrix,2);off_matrix(off_zero,:) = [];
fit_variance_explained = cell(36,1); % variance of the rSC's response explained by the response of 36 ROIs in the left-hemisphere
for r = 1:36    
    on_test_orig = test_on{r}(:,6:25);
    off_test_orig = test_off{r}(:,6:25);    
    on_test_orig(isnan(on_test_orig))=0; on_test_orig(on_zero,:) = [];
    off_test_orig(isnan(off_test_orig))=0; off_test_orig(off_zero,:) = [];    
    % 1. split trials into 10 portions; 9 portions for fitting; 1 portion for testing:
    rng('shuffle');
    N = size(on_matrix,1);
    F = size(off_matrix,1);
    on_trial_id = randperm(size(on_matrix,1));
    off_trial_id = randperm(size(off_matrix,1));    
    fit_variance_explained{r} = zeros(10,1);
    for p=1:10
        if p<10
            on_tid = on_trial_id(1+(p-1)*round(N/10):(p-1)*round(N/10)+round(N/10));
            off_tid = off_trial_id(1+(p-1)*round(F/10):(p-1)*round(F/10)+round(F/10));
        else
            on_tid = on_trial_id(1+(p-1)*round(N/10):end);
            off_tid = off_trial_id(1+(p-1)*round(F/10):end);
        end
        on_tarpredict = on_matrix(on_tid,:);
        off_tarpredict = off_matrix(off_tid,:);
        on_target = on_matrix(setdiff(on_trial_id,on_tid),:);
        off_target = off_matrix(setdiff(off_trial_id,off_tid),:);
        
        on_testpredict=on_test_orig(on_tid,:);
        off_testpredict=off_test_orig(off_tid,:);
        on_test = on_test_orig(setdiff(on_trial_id,on_tid),:);
        off_test = off_test_orig(setdiff(off_trial_id,off_tid),:);
        
        X_ = zeros(41, 41);
        for j = 1:20
            X_(j,1) = sum(on_test(:,j));
            X_(j,j+1) = size(on_test,1);
        end
        for j=21:40
            X_(j,1) = sum(off_test(:,j-20));
            X_(j,j+1) = size(off_test,1);
        end
        X_(41,1) = sum(sum(on_test.^2))+sum(sum(off_test.^2));
        for i = 1:20
            X_(41,i+1) = sum(on_test(:,i));
        end
        for i= 21:40
            X_(41,i+1) = sum(off_test(:,i-20));
        end
        
        Y_ = zeros(41,1);
        for i = 1:20
            Y_(i) = sum(on_target(:,i));
        end
        for i = 21:40
            Y_(i) = sum(off_target(:,i-20));
        end
        Y_(41) = sum(sum(on_test.*on_target))+sum(sum(off_test.*off_target));
        
        V_ = X_\Y_;
        
        fit_on_SC = repmat(V_(2:21)',size(on_tarpredict,1),1)+V_(1)*on_testpredict;
        fit_off_SC = repmat(V_(22:41)',size(off_tarpredict,1),1) + V_(1)*off_testpredict;
        
        TSS = sum(sum((off_tarpredict-mean([off_tarpredict;on_tarpredict])).^2)) + sum(sum((on_tarpredict-mean([off_tarpredict;on_tarpredict])).^2));
        RSS = sum(sum((off_tarpredict-fit_off_SC).^2))+sum(sum((on_tarpredict-fit_on_SC).^2));
        fit_variance_explained{r}(p) = 1-RSS/TSS;
    end
end

%% 3. Simultaneous fiber photometry and fMRI
% Note: function "coherencyc" is from the Chronux: http://chronux.org/
clear;
path = 'YOURPATH/fiber photometry fMRI';
animals = dir(path);
fMRI_30 = [];
Ca_30 = [];
fMRI_20 = [];
Ca_20 = [];
for i = 1:size(animals,1)    
    id = strfind(animals(i).name, 'C');
    if ~isempty(id)
        animal_id = animals(i).name(id+1:end); 
        cd(fullfile(path,['SC',num2str(animal_id)]));
        animal_folder = cd;        
        scan_dates = dir(cd);
        for d = 1:length(scan_dates)-2
            epi_folder = fullfile(animal_folder, scan_dates(d+2).name); 
            cd(epi_folder);
            EPIs = dir(cd);
            for e = 1:length(EPIs)-2
                cd(fullfile(epi_folder, EPIs(e+2).name));%
                load('trigger'); load('fiber_mask'); load('BOLD'); load('GCaMP');
                if epi_frames_off(2)-epi_frames_off(1)>45                    
                    fMRI_30 = [fMRI_30, mean(BOLD(fiber_mask,11:600))'];
                    Ca_30 = [Ca_30, GCaMP(11:600)];
                else
                    fMRI_20 = [fMRI_20, mean(BOLD(fiber_mask,11:600))'];
                    Ca_20 = [Ca_20, GCaMP(11:600)];
                end 
            end            
        end
    end    
end            

% calculate the magnitude-squared coherence
params.pad = 0;
params.tapers = [3,5];
params.err = [2, 0.05];
params.trialave = 1;

[C_30,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(Ca_30,fMRI_30,params);
 
[C_20,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(Ca_20,fMRI_20,params);










