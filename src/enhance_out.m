function [ I_t I_r ] = enhance_out(I_t, I_r, padding, I_in)

% Post-processing to enhance the results
% Subtract min value to center the image
I_t_retrieve = I_t(padding + 1: end - padding - 1, padding + 1: end - padding);
I_t = I_t - min(I_t_retrieve(:));

I_r_retrieve = I_r(padding + 1: end - padding - 1, padding + 1: end - padding);
I_r = I_r - min(I_r_retrieve(:));

% Match the global information of the original image with that of the
% transmitted image
sig = sqrt( sum( (I_in - mean(I_in) ).^2 ) / sum( (I_t-mean(I_t(:)) ).^2 ) );
I_t = sig * ( I_t - mean(I_t(:)) ) + mean( I_in(:) );
