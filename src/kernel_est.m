function [dx dy c ] =  kernel_est(I) 

    I = rgb2gray(I);
    auto_corr_map = auto_corr_on_laplacian(I);
    f = @(x) max(x(:));
    local_maxima = nlfilter(auto_corr_map,[5 5],f);
    f2 = @(x) max(x(x~=max(x)));
    second_local_maxima = nlfilter(auto_corr_map,[5 5],f2);

    [size_x, size_y] = size(local_maxima);

    local_maxima((size_x+1)/2 - 4: (size_x+1)/2 + 4, (size_y+1)/2 - 4: (size_y+1)/2 + 4) = 0;
    second_local_maxima((size_x+1)/2 - 4: (size_x+1)/2 + 4, (size_y+1)/2 - 4: (size_y+1)/2 + 4) = 0;

    difference = local_maxima - second_local_maxima;
    threshold = 70;
    indices = find(difference < threshold);
    local_maxima(indices) = 0;
    [global_maxima, dk] = max(local_maxima(:));
    [dk_y, dk_x] = ind2sub(size(local_maxima),dk);
    dk_y = floor((size_x)/2 + 1 - dk_y)
    dk_x = floor((size_y)/2 + 1 - dk_x) 

    ck = estimate_attenuation(I, dk_x, dk_y)

    % Visualization
    figure;
    subplot(1, 2, 1);
    imagesc(local_maxima);
    colorbar;
    subplot(1, 2, 2);
    imagesc(second_local_maxima);
end