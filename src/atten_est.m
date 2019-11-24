function c = est_attenuation(I, dx, dy)
  % Function estimates attenuation factor c_k using
  % Harris corner detector

  fprintf('Estimating attenuation factor...\n');

  % Harris Corner Detector
  cns = corner(I);
  patch_size = 5;
  %m = [];

  attn = zeros(size(cns,1),1);
  weight = zeros(size(cns,1),1);

  for i = 1 : size(cns, 1)

    % Get a patch and a patch at spatial offset d_k
    patch_1 = get_patch(I, cns(i, 1), cns(i, 2), patch_size);
    patch_2 = get_patch(I, cns(i, 1) + dx, cns(i, 2) + dy, patch_size);

    if isempty(patch_1) || isempty(patch_2)
      continue;
    end

    % m = [m; patch_1 patch_2;];

    % Calculating variance of original patch and spatial offset patch
    patch_1 = patch_1(:);
    patch_2 = patch_2(:);
    % Normalizing patch wrt to mean
    patch_1 = patch_1 - mean(patch_1);
    patch_2 = patch_2 - mean(patch_2);
  
    attn(i) = (max(patch_2) - min(patch_2)) / (max(patch_1) - min(patch_1));

    if (attn(i) > 0) && (attn(i) < 1)
      score = sum(patch_1.*patch_2) / sqrt(sum(patch_1.^2)) * sqrt(sum(patch_2.^2)) ;
      weight(i) = exp(-score / (2*0.2^2));
    end

  end

  % Apply normalization factor on the formula for c
  c = sum(weight.*attn) / sum(weight);

end

function patch = get_patch(I, x, y, patch_size)
  % Using a 5x5 patch
  if (x > patch_size) && (y > patch_size) && (x < size(I, 2) - patch_size) && (y < size(I, 1) - patch_size)
    patch = I(y - patch_size: y + patch_size, x - patch_size: x + patch_size);
  else
    patch = [];
  end
end
