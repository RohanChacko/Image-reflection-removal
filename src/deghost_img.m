function deghost_image(img_path)

% Dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch-GMM prior
addpath('epllcode');
% Bounded-LBFGS Optimization
addpath('lbfgsb/lbfgsb3.0_mex1.2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_in = im2double(imread(img_path));
fprintf('Size of image: %d %d %d\n', size(I_in));

% Ghosting kernel consists of spatial shift vector d_k
% and attenuation factor c_k
fprintf('Estimating Ghosting kernel...\n');

% Estimate Ghosting kernel from input image
[configs.dx configs.dy configs.c] = kernel_est(I_in);
fprintf('Done.\n');

% Padding size : needs to be larger than the spatial shift.
configs.padding = ceil(norm([configs.dx configs.dy]))+10;

[configs.h configs.w channels] = size(I_in);
configs.num_px = configs.h * configs.w;

% Patch-GMM Optimization

for i = 1 : channels
  fprintf('Applying Patch GMM on channel %d...\n', i);
  configs.ch=i;

  % Applying Patch-GMM Optimization to each channel separately
  [I_t_k I_r_k ] = patch_gmm(I_in(:,:,i), configs);

  % Post-processing to enhance the results
  [I_t(:,:,i), I_r(:,:,i)] = enhance_out(I_t_k, I_r_k, configs.padding, I_in(:,:,i));

end
fprintf('Done.\n');

%% Output results
% <filename>_input.png        -> original input image
% <filename>_transmission.png -> transmission layer
% <filename>_reflection.png   -> reflection layer
fprintf('Saving results\n');
tp = split(img_path, '/');
tp = split(tp(end), '.');

imwrite(I_in, char(strcat(tp(1),'_input.png')));
imwrite(I_t, char(strcat(tp(1),'_transmission.png')));
imwrite(I_r, char(strcat(tp(1),'_reflection.png')));
