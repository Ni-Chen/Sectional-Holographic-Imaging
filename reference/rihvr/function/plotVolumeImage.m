function plotVolumeImage(intensity)

volume_size = size(intensity);
num_planes = volume_size(3);
image_size = volume_size(1:2);
plot_size = image_size;

N = image_size(1);
ySchnitt = 0;  % position of xz-cut out (counted from the center y)
zSchnitt = -25;   % position of xy-cut out (counted from the center z)
T1 = maketform('affine', [1 0; 0 1; 0 0]);
inverseFcn = @(X, t) [X repmat(t.tdata, [size(X, 1) 1])];
T2 = maketform('custom', 3, 2, [], inverseFcn, ceil(N/2) + ySchnitt);
Tc = maketform('composite', T1, T2);
R3 = makeresampler({'cubic', 'nearest', 'nearest'}, 'fill');

% close all;
plotPicture_flipped(intensity, plot_size(1), plot_size(2), Tc, R3, 0, '3D', 0.001);
% title('Inverse');
drawnow();

fig = gcf;
fig.Color = 'white';
axis normal;
xticks([]);
yticks([]);
zticks([]);
xlabel('');
ylabel('');
zlabel('');
colorbar('off');
axis equal;
frame = getframe(fig);
im = frame2im(frame);
imwrite(im, 'Outputs/volume_inverse.png');
% close all;