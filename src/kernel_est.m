function [dx dy c] =  kernel_est(I_in)
  % Estimate kernel parameters :
  % 1. d_k : spatial shift vector
  % 2. c_k : attenuation factor

  fprintf('Estimating spatial shift offset...\n');

  % Converting image to grayscale
  I_in = rgb2gray(I_in);

  OFFSET = 50;
  laplacian_filter = [0 -1 0; -1  4 -1; 0 -1 0];

  % Applying laplacian on the image
  lp_output = imfilter(I_in, laplacian_filter);
  % Generaing autocorrelation map
  auto_corr = xcorr2(lp_output);
  [corr_x , corr_y] = size(auto_corr);
  auto_corr = auto_corr(floor((corr_x+1)/2)-OFFSET:floor((corr_x+1)/2)+OFFSET, floor((corr_y+1)/2)-OFFSET:floor((corr_y+1)/2)+OFFSET);

  % Extracting the first and second local maxima in each 5x5 neighborhood
  f1 = @(x) max(x(:));
  f_loc_max = nlfilter(auto_corr, [5 5], f1);
  f2 = @(x) max(x(x~=max(x)));
  s_loc_max = nlfilter(auto_corr, [5 5], f2);

  % Removing local maxima within 4 pixels of origin for robust estimation
  [size_x, size_y] = size(f_loc_max);
  f_loc_max((size_x+1)/2 - 4: (size_x+1)/2 + 4, (size_y+1)/2 - 4: (size_y+1)/2 + 4) = 0;
  s_loc_max((size_x+1)/2 - 4: (size_x+1)/2 + 4, (size_y+1)/2 - 4: (size_y+1)/2 + 4) = 0;

  % Discard local maxima in neighborhoods where the first & second maxima
  % are closer than a predefined threshold to remove local maxima
  % caused due to locally flat or repetitive structures. Here threshold : 70

  threshold = 70;
  %thresh_diff = f_loc_max - s_loc_max;
  %indx = find(thresh_diff < threshold);
  %f_loc_max(indx) = 0;
  f_local_max(f_loc_max - s_loc_max < 70) = 0;

  % Select the largest maxima from the local maximas as Ghosting distance
  [global_maxima, dk] = max(f_loc_max(:));
  [dk_y, dk_x] = ind2sub(size(f_loc_max), dk);

  % Centering the result wrt to origin
  dy = floor((size_x)/2 + 1 - dk_y)
  dx = floor((size_y)/2 + 1 - dk_x)

  c = atten_est(I_in, dx, dy);
end
