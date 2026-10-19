function [k_min, N_j, info] = recubrimiento_Nj(s, p_max, c)
% From the geometry of the covering alone (conditions on n and
% Nj, which depend neither on f nor on the particular value of phi), it returns:
% k_min : the first k such that N=(2^k+1)^s satisfies condition (1):
% N/N_j>=n_0^s.
% N_j: smallest (nominal) number of nodes per patch with N_j >= max{Nj2, Nj3}.
% INPUTS:
% s: dimension.
% p_max: largest odd p used, 5 by default. It only appears in (3) through
% L=c(s+(p_max-1)/2,s)).
% c: slack in (2) and (3), 1.75 by default. Lowering it reduces N_j.
% OUTPUT: info struct with n0, lambda_s, alpha0, L, Nj2, Nj3, c, p_max. To
% extract one of those variables from the struct we have to write, for example,
% "info.n0".
    if nargin < 2 || isempty(p_max), p_max = 5;   end % if we pass it
      % fewer than two input arguments, or the second input is an empty
      % vector, it sets p_max=5 by default.
      if nargin < 3 || isempty(c), c = 7/4; end %if we pass it fewer than three
      % arguments, or the third one is an empty vector, it sets c=7/4.
      lambda_s = pi^(s/2) / gamma(s/2 + 1);% volume of the unit ball
      n0 = floor(1/(1 - sqrt(s)/2)) + 1; % condition on n.
      alpha0 = (n0 - 1)/n0;
      m = (p_max + 1)/2;
      L = nchoosek(s + (m - 1), s); % dim Pi_{m-1}^s, for the largest p.

      base = c^s / (lambda_s*alpha0^s);
      Nj2 = base;
      Nj3 = L * base;
      N_j = ceil(max(Nj2, Nj3));

      k_min = 1;
      while round( (2^k_min + 1) / N_j^(1/s) ) < n0
          k_min = k_min + 1;
      end

      info = struct('n0', n0, 'lambda_s', lambda_s, 'alpha0', alpha0, ...
                    'L', L, 'Nj2', Nj2, 'Nj3', Nj3, 'c', c, 'p_max', p_max);
end
