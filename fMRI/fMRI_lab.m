%%
addpath('matlab\');

%%
 

fMRI = MRIread('conv.nii');
t1 = MRIread('t1.nii');
% register volume i to 1st volume
fMRIregistered = fMRI;
for i=1:size(fMRI,4)
    [~,fMRIregistered.vol(:,:,:,i) ] = ...
        imregmoment(fMRIregistered.vol(:,:,:,i),fMRIregistered.vol(:,:,:,1));
end

x1len = size(fMRIregistered.vol, 1);
x2len = size(fMRIregistered.vol, 2);
x3len = size(fMRIregistered.vol, 3);
volsize = x1len*x2len*x3len;
x4len = size(fMRIregistered.vol, 4);

%%

figure;
temp = permute(fMRIregistered.vol, [1,2,4,3]);
imshow3D(temp(:,:,:,1));
clear temp;

%%
fMRIvr = reshape(fMRIregistered.vol, [volsize,x4len]);
fMRIvr_mean = mean(fMRIvr,2);
fMRIvr_nlized = fMRIvr - fMRIvr_mean;

fMRIvr_std = std(fMRIvr_nlized,0,2);

signal_ratio = fMRIvr_std ./ fMRIvr_mean;
signal_ratio(isnan(signal_ratio)) = 0;

disp(mean(signal_ratio((round(volsize/2)-20):(round(volsize/2)+20))));
