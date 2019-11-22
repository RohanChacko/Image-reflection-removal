function [dx dy c ] =  kernel_est(I_in) 

  I_in = rgb2gray(I_in);
  Laplacian=[0 -1 0; -1 4 -1; 0 -1 0];
  resp = imfilter(I_in, Laplacian);
  auto_corr = xcorr2(resp, resp);
  bdry = 370; 
  auto_corr = auto_corr(bdry:end-bdry, bdry:end-bdry);
  
  max_1 = ordfilt2(auto_corr, 25, true(5));
  max_2 = ordfilt2(auto_corr, 24, true(5));

  auto_corr(end/2 - 4 : end/2 + 4, end/2 - 4 : end/2+4)=0;
  candidates = find((auto_corr == max_1) & ((max_1 - max_2)>70));
  candidates_val = auto_corr(candidates);

  cur_max = 0;
  dx = 0; 
  dy = 0;
  offset = size(auto_corr)/2 + 1;
  for i = 1 : length(candidates)
    if (candidates_val(i) > cur_max)  
      [dy dx] = ind2sub(size(auto_corr), candidates(i)); 
      dy = dy - offset(1);
      dx = dx - offset(2);
    end
  end
  c = est_attenuation(I_in, dx, dy);
