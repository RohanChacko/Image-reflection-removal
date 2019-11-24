function [I_t I_r ] = patch_gmm(I_in, configs)

  % Setup for patch-based reconstruction
  h = configs.h;
  w = configs.w;
  c = configs.c;
  dx = configs.dx;
  dy = configs.dy;
  lambda = 1e6;

  % patch size for patch-GMM
  psize = 8;

  % Parameters for half-quadratic regularization method
  configs.beta_factor = 2;
  configs.beta_i = 200;
  configs.dims = [h w];

  niter = 5;
  beta  = configs.beta_i;

  % Identity matrix
  Id_mat = speye(h*w, h*w);

  % Ghosting kernel K
  k_mat = get_k(h, w, dx, dy, c);

  % Operator that maps an image to its ghosted version
  A = [Id_mat k_mat];

  num_patches = (h - psize + 1) * (w - psize + 1);

  mask = merge_two_patches(ones(psize^2, num_patches),...
              ones(psize^2, num_patches), h, w, psize);

  % Use non-negative constraint
  configs.non_negative = true;

  % Setup for GMM prior
  load GSModel_8x8_200_2M_noDC_zeromean.mat
  excludeList=[];

  % Initialize
  [I_t_i I_r_i ] = grad_irls(I_in, configs);

  % faster option, but results are not as good.
  % [I_t_i I_r_i ] = grad_lasso(I_in, configs);


  % Apply patch gmm with the initial result.
  % Create patches from the two layers.
  est_t = im2patches(I_t_i, psize);
  est_r = im2patches(I_r_i, psize);


  for i = 1 : niter
    fprintf('Optimization iteration : %d...\n', i);

    % Merge the patches with bounded least squares
    sum_piT_zi = merge_two_patches(est_t, est_r, h, w, psize);
    z = lambda * transpose(A) * I_in(:) + beta * sum_piT_zi;

    % Non-negative optimization by L-BFGSB
    opts = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 50);
    opts.printEvery = 50;

    sum_zi_2 = norm(est_t(:))^2 + norm(est_r(:))^2;
    fcn = @(x)( lambda * norm(A*x - I_in(:))^2 + ...
        beta*( sum(x.*mask.*x - 2 * x.* sum_piT_zi(:)) + sum_zi_2));

    f_handle = @(x)(lambda * transpose(A) * (A*x) + beta * (mask.*x));
    grad = @(x)(2*(f_handle(x) - z));
    fun = @(x)fminunc_wrapper( x, fcn, grad);

    l = zeros(numel(z), 1);
    u = ones(numel(z), 1);

    [out, ~, info] = lbfgsb(fun, l, u, opts );

    out = reshape(out, h, w, 2);
    I_t = out(:,:,1);
    I_r = out(:,:,2);

    % Extract est_t and est_r by restoring patches using the patch based prior
    est_t = im2patches(I_t, psize);
    est_r = im2patches(I_r, psize);

    noiseSD=(1/beta)^0.5;
    [est_t t_cost]= aprxMAPGMM(est_t, psize, noiseSD, [h w], GS, excludeList);
    [est_r r_cost]= aprxMAPGMM(est_r, psize, noiseSD, [h w], GS, excludeList);

    beta = beta*configs.beta_factor;

  end
