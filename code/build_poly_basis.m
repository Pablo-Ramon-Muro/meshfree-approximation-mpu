function P = build_poly_basis(pts, p, center, rho)
    % Builds the polynomial basis Pi_{m-1}^s evaluated at the points 'pts',
    % which is used to build the local interpolant with phi=r^p. P is a
    % matrix whose first column is the first basis vector evaluated at each
    % point of pts, the second column the second vector evaluated at all
    % the pts, and so on.
    % s = 2 (dimension of the space: x, y)
    % m = (p+1)/2. Therefore:
    % If p=1 -> m=1 -> Polynomials of degree 0 (constant)
    % If p=3 -> m=2 -> Polynomials of degree <= 1 (linear)
    % If p=5 -> m=3 -> Polynomials of degree <= 2 (quadratic)
    % If p=0 -> m=0 -> Gaussian (does not need polynomials)
    %
    % The basis is recentred and rescaled to the patch: xt=(x-cx)/rho,
    % yt=(y-cy)/rho. It spans the same space Pi_{m-1}^s.
    x = (pts(:, 1) - center(1)) / rho;
    y = (pts(:, 2) - center(2)) / rho;
    N_pts = length(x);

      if p == 0
          P = [];                                                  % Gaussian
      elseif p == 1
          P = ones(N_pts, 1);                                      % L = 1
      elseif p == 3
          P = [ones(N_pts, 1), x, y];                          % L = 3
      elseif p == 5
          P = [ones(N_pts, 1), x, y, x.^2, x.*y, y.^2];        % L = 6
      else
          error('Valor de p no soportado. Elige p = 1, 3, 5 o 0 (Gauss).');
      end
end
