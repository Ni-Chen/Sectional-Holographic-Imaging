close all; clear; clc;

load('beads_TwIST_result');

bp_reobj1 = reobj_raw(:,:,1);
bp_reobj2 = reobj_raw(:,:,2);
temp1 = abs(bp_reobj1);
temp2 = abs(bp_reobj2);
% bp_reobj1 = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
% bp_reobj2 = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
bp_reobj1 = mat2gray(temp1, [0 255]);
bp_reobj2 = mat2gray(temp2, [0 255]);

figure; imshow(abs(bp_reobj1),[],'Border','tight', 'InitialMagnification', 100);  colormap(hot);%#    a border at full magnification
set(gca,'LooseInset',get(gca,'TightInset')); print('bp_reobj1.png','-dpng','-r300');
figure; imshow(abs(bp_reobj2),[],'Border','tight', 'InitialMagnification', 100);  colormap(hot);%#    a border at full magnification
set(gca,'LooseInset',get(gca,'TightInset')); print('bp_reobj2.png','-dpng','-r300');

twist_reobj1 = reobj_deconv(:,:,1);
twist_reobj2 = reobj_deconv(:,:,2);
temp1 = abs(twist_reobj1);
temp2 = abs(twist_reobj2);
% twist_reobj1 = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
% twist_reobj2 = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
twist_reobj1 = mat2gray(temp1, [0 255]);
twist_reobj2 = mat2gray(temp2, [0 255]);

figure; imshow(abs(twist_reobj1),[],'Border','tight', 'InitialMagnification', 100);  colormap(hot);%#    a border at full magnification
set(gca,'LooseInset',get(gca,'TightInset')); print('twist_reobj1.png','-dpng','-r300');
figure; imshow(abs(twist_reobj2),[],'Border','tight', 'InitialMagnification', 100);  colormap(hot);%#    a border at full magnification
set(gca,'LooseInset',get(gca,'TightInset')); print('twist_reobj2.png','-dpng','-r300');

load('blind_tikreg_recon');
tikreg_reobj1 = imresize(reobj2, [512 512]);
tikreg_reobj2 = imresize(reobj1, [512 512]);
temp1 = abs(tikreg_reobj1);
temp2 = abs(tikreg_reobj2);

tikreg_reobj1 = mat2gray(temp1, [0 255]);
tikreg_reobj2 = mat2gray(temp2, [0 255]);

figure; imshow(tikreg_reobj1,[],'Border','tight', 'InitialMagnification', 100);  colormap(hot);%#    a border at full magnification
set(gca,'LooseInset',get(gca,'TightInset')); print('tikreg_reobj1.png','-dpng','-r300');
figure; imshow(tikreg_reobj2,[],'Border','tight', 'InitialMagnification', 100);  colormap(hot);%#    a border at full magnification
set(gca,'LooseInset',get(gca,'TightInset')); print('tikreg_reobj2.png','-dpng','-r300');



