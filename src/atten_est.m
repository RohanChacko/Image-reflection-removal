function c = est_attenuation(I, dx, dy)
  % Function estimates attenuation factor c_k using
  % Harris corner detector

  % Extracting relevant corners
  cns = corner(I);
  hw = 5;
  m = [];

  score=zeros(size(cns,1),1);
  attn=score*0;
  w=score*0;

  for i = 1 : size(cns,1)

    % Get a patch and a patch at spatial offset d_k
    p1 = get_patch(I, cns(i,1), cns(i,2), hw);
    p2 = get_patch(I, cns(i,1) + dx, cns(i,2) + dy, hw);

    if ~isempty(p1) && ~isempty(p2)
      m = [m; p1 p2;];

      % Calculate variances of patch and patch+d_k
      p1 = p1(:);
      p1 = p1 - mean(p1);
      p2 = p2(:);
      p2 = p2 - mean(p2);

      score(i) = sum(p1.*p2) / sum(p1.^2).^0.5 / sum(p2.^2).^0.5;
      attn(i) = (max(p2) - min(p2)) / (max(p1) - min(p1));

      if (attn(i) <1) && ( attn(i) > 0)

        w(i) = exp(-score(i)/(2*0.2^2));
      end
    end
  end

  % Apply normalization factor on the formula for c
  c = sum(w.*attn)/sum(w);

end

function p = get_patch(I, x, y, hw)
  if (x > hw) && (y > hw) && (x < size(I, 2) - hw) && (y < size(I, 1) - hw)
    p= I(y - hw: y + hw, x - hw: x + hw);
  else
    p=[];
  end
end