function show_brain_map_highres(map,map_wholefov,min,max,cluster_size,selected_slices,option)

ana=spm_read_vols(spm_vol(which('standard_anatomy_t2.nii')));
ana = ana(:,:,selected_slices);
ana(3:256,:,:) = ana(1:254,:,:);
brain_mask=spm_read_vols(spm_vol(which('brain_mask_64x64.nii')));
% wm_csf_mask = spm_read_vols(spm_vol(which('WM_CSF_mask_64x64.nii')));
brain_mask = brain_mask(:,:,selected_slices);
% wm_csf_mask = wm_csf_mask(:,:,selected_slices)>0;
% define_colormap;

if isempty(map_wholefov)
    map3d = zeros(64,64,length(selected_slices));
    map3d(brain_mask>0) = map(:);
    map = map3d;
    % map(wm_csf_mask) = nan;
else
    map = reshape(map_wholefov,64,64,[]);
end

ana = ana(61:200,256:-1:129,:);

for i_slice = 1:size(ana,3)
    x = ceil(i_slice/5);
    y = rem(i_slice,5);
    if y==0
        y=5;
    end
    ana_flat(128*(x-1)+1:128*x,140*(y-1)+1:140*y) = ana(:,:,i_slice)';
end
ana_flat(end+30,:,:) = 0;

map = imresize(map, 4, 'nearest');
map = map(61:200,256:-1:129,:);
for i = 1:size(map,3)
    x = ceil(i/5);
    y = rem(i,5);
    if y==0
        y=5;
    end
    map_flat(128*(x-1)+1:128*x,140*(y-1)+1:140*y) = map(:,:,i)';
end
map_flat(end+30,:,:) = 0;

if option == 1 % symmetric color range
    mask = abs(map_flat)>min;
    mask = bwareaopen(mask,cluster_size,6);
    map_flat(~mask) = nan;
%     
    h = figure;
    h.Position = [0,0,800,800];
    rat_fmri_imoverlay(ana_flat,map_flat,...
        [-max max],[0,5800],'jet',0.8,h);
    colormap('jet');
else
    mask = map_flat>min;
    mask = bwareaopen(mask,cluster_size,8);
    map_flat(~mask) = nan;
%     
    h = figure;
    h.Position = [0,0,800,800];
    rat_fmri_imoverlay(ana_flat,map_flat,...
        [min max],[0,5800],'parula',0.8,h);
end
% colormap(map_b);

bregmas = -9.2:1:4.8;

text(5,36,['Bregma ',num2str(-10.2),' mm'],'FontSize',12,'FontWeight','bold','Color','w');

ys = 48 + (1:140:620);
xs = 140 + (1:128:380);
% xs(end) = 128;


for i = 1:length(bregmas)
    x = floor((i-1)/5)+1;
    y = rem(i,5);
    if y == 0
        y = 5;
    end
    text(ys(y),xs(x),[num2str(bregmas(i)),' mm'],'FontSize',12,'FontWeight','bold','Color','w');
end

end
