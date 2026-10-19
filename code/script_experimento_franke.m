% =========================================================================
% SCRIPT: MPU for phi=r^p, with p = 1, 3, 5).
% Adapted for different weight functions (Wendland).
% =========================================================================
%   - Sweep in N=(2^k+1)^s with N_j fixed and M=M(N)
%   - eps of the Gaussian chosen automatically (worst patch of the largest N,
% condition 4).
%   - phi_rho= (r/rho)^p (Lemma 1.60) .
%   - polynomial basis rescaled on each patch.
%   - one convergence table per combination (with h and min(RCOND)).
%
% All the geometry that depends only on N (tree of
% nodes, DM_eval, indices of nodes per patch and the weight matrices SEM)
% is computed only once per N, outside the loop over the combinations.
% Besides, to make it as efficient as possible, the local interpolant is
% computed once per Phi and reused for the weights
% (it does not depend on the weight).
% =========================================================================

clear; clc;

% --- 1. PARAMETER SETUP ---
% n and wep follow from N and N_j.
% 'ep' (Gaussian shape) is computed automatically (section 2c).
neval = 100;                  % Evaluation points per axis (total: neval^2)
rcond_target = 1e-13;        % conditioning threshold (RCOND)

% We define the test function: Franke's in this case.
f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
func_test_franke = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);

% --- 2. GENERATION OF POINTS ---
% 2a. Evaluation points (regular uniform grid)
% neval is fixed => epoints and exact do not depend on N: computed only once.
grid_eval = linspace(0, 1, neval);
[xe, ye] = meshgrid(grid_eval);
epoints = [xe(:), ye(:)]; %neval^2x2 matrix where each row is an evaluation
% point.
exact = func_test_franke(epoints(:,1), epoints(:,2)); % Exact solution:
% Franke's function evaluated at the evaluation points.

% --- 2b. CHOICE OF N_j AND RANGE OF THE SWEEP ---
% N_j and k_min come out of the function recubrimiento_Nj, which
% applies the conditions on the covering (they depend only on s, on the
% largest p used and on the constant c=1.75).
s = 2;     % dimension (Franke on [0,1]^2)
[k_min, N_j, recinfo] = recubrimiento_Nj(s, 5, 1.75); % p_max = 5, c = 1.75
% and we get k_min, which will give us N_min. It also gives us the fixed N_j we
% are after, and recinfo collects n0, L, Nj2, Nj3, c.
fprintf(['N_j (recubrimiento) = %d    (n0=%d, L=%d, Nj2=%.2f, Nj3=%.2f,' ...
  ' c=%.2f)\n'], N_j, recinfo.n0, recinfo.L, recinfo.Nj2, recinfo.Nj3, ...
  recinfo.c);%gives us the value of N
k_max = 10;    % N up to (2^10+1)^s
k_list = k_min:k_max;
nK = numel(k_list); %length of k_list
fprintf('Barrido: k=%d..%d -> N=%d..%d\n', k_min, k_max, ...
  (2^k_min+1)^s, (2^k_max+1)^s);

% --- 2c. AUTOMATIC SELECTION OF eps_gauss (real worst patch for N_max,
% the largest N)
% The Gaussian conditioning is dominated by the closest pair of nodes
% inside a patch. We locate the most tightly packed nodes for N_max,
% we look at the patches containing them, we take the one with the smallest rcond and
% we tune eps by bisection so that rcond(worst) >= rcond_target.
% rcond grows with eps. For smaller N the conditioning only improves.
ep_native = 9/sqrt(2); % epsilon threshold
N_max   = (2^k_max + 1)^s;
n_max = round((N_max/N_j)^(1/s));
if n_max < recinfo.n0 %checks that N_max is large enough
  % to meet the bounds required by the covering.
    error('MPU:sinCobertura', ...
        ['Ni N_max=%d alcanza cobertura: n_max=%d < n0=%d ' ...
          '(alpha=(n-1)/n <= sqrt(s)/2). Sube k_max o baja N_j.'], ...
        N_max, n_max, recinfo.n0);
end
rho_max = 1/n_max;
ds_max = haltonseq(N_max, s);%Halton nodes for N_max
tree_mx = createns(ds_max, 'NSMethod', 'kdtree');
gpu_mx = linspace(0, 1, n_max);
[~, Dnn] = knnsearch(ds_max, ds_max, 'K', 2);%in order not to compute every
% distance, it looks for the two nearest neighbours (Halton nodes) of each
% Halton node, and Dnn stores in its first column the nearest
% distances (0, since it is the node with itself) and in the second column the
% smallest distance to the nearest node (which node that is,
% that is, its index, is not stored). Each row of Dnn gives the two
% nearest nodes to a node (following an order).
[~, ord] = sort(Dnn(:,2), 'ascend');%it keeps the 2nd column of Dnn, sorts
% the distances from smallest to largest and keeps the index of the nodes
% ordered from smallest to largest proximity to another node.
cand_nod = ds_max(ord(1:min(300,end)), :);%we keep the 300 nodes
% with the smallest distances (we store their coordinates in an N_max x 2 matrix).
gi         = @(v) min(max(round(v*(n_max-1)) + 1, 1), n_max);
cand_ctr = unique([gpu_mx(gi(cand_nod(:,1))).', ...
  gpu_mx(gi(cand_nod(:,2))).'], 'rows');%what the function gi does is take the
% interval [0,1] and, by multiplying v*(n_max-1), rescale it to [0,n_max-1),
% round it to the nearest n and add one, so as to take an element of
% {1,...,n_max}. This will be the index of the centre
% nearest to a node in that dimension (not the total distance, but in that
% dimension). We apply that coordinate by coordinate to the 300 worst
% nodes and we get two columns. Each column corresponds to the element of
% {1,...,n_max} on that axis. Finally gpu_max(coordinate) tells us which
% x coordinate (y respectively) of a centre that index corresponds to.
% unique by rows removes the repeated rows. So cand_ctr gives us the
% coordinates of the patches that will be the worst conditioned.
worst_rc = Inf; dist_rep = []; nrep = 0; %we initialise the worst rcond at a
% large number (it will come down as we iterate), as well as the distance
% matrix associated with the patch and the number of nodes of the patch.
for cc = 1:size(cand_ctr,1)
     [ic,~] = rangesearch(tree_mx, cand_ctr(cc,:), rho_max);%it stores the
     % index of the nodes falling inside the cc-th problematic patch.
     nod = ds_max(ic{1}, :);%nod stores the x,y coordinates of the nodes
     %inside the cc-th patch.
     if size(nod,1) < 2, continue; end %if there is only one node in the
     % cc-th patch, what comes next is ignored and the for loop starts again
     % on the next patch.
     Dloc = DistanceMatrix(nod, nod);%it computes the distance matrix on that
     % patch.
     rcref = rcond(exp(-(ep_native*Dloc).^2));%it gives us the rcond of
     % Phi(distmatrix).
     if rcref < worst_rc, worst_rc = rcref; %if the new conditioning is
       % worse than the previous one, we store this new rcond (as the worst),
       % along with its distance matrix and the number of nodes in the patch.
       dist_rep = Dloc; nrep = size(nod,1);
     end
end
gr = @(ep) rcond(exp(-(ep*dist_rep).^2));
if gr(ep_native) >= rcond_target
     ep = ep_native;%if, for the initial epsilon, the worst conditioned
     %local system is as well conditioned as or better than rcond=10^(-13), we
     % keep that epsilon. Then initially the length of the interval
     % [ep_lo,ep_hi] is about 10^7. If we halve it 80 times,
     % that is of the order of 10^(-24), which multiplied by 10^7 leaves us
     % a window of length 10^(-17), while the precision in MATLAB is
     % 10^(-16), so it is enough.
else
    ep_lo = ep_native; ep_hi = 2*ep_native;
    while gr(ep_hi) < rcond_target, ep_hi = 2*ep_hi; if ep_hi>1e7, break;
    end, end %it gives us a floor, ep_lo (the epsilon whose rcond is smaller than
  % 10^(-13)) and a ceiling, ep_hi (this one has rcond>=10^(-13)).
    for it = 1:80
        em = 0.5*(ep_lo+ep_hi);
        if gr(em) >= rcond_target, ep_hi = em; else, ep_lo = em; end
    end%bisection is carried out in such a way that 80 iterations are enough
        % because, in the extreme case ep_hi=10^7, initially
        % the length of the interval [ep_lo,ep_hi] is
        % about 10^7, and if we halve it 80 times,
    % that is of the order of 10^(-24), which multiplied by 10^7 leaves us
    % a window of length 10^(-17), while the precision in MATLAB is
    % 10^(-16), so it is enough.
    ep = ep_hi;%now I do keep ep_hi, as the smallest ep for
    % which it will be well conditioned.
end
ep = max(ep, ep_native); % safeguard so that it is not smaller than the initial
% epsilon.
fprintf(['eps_gauss = %.4f (peor parche a N_max: %d nodos,' ...
  ' rcond=%.2e)\n\n'], ep, nrep, gr(ep));
clear ds_max tree_mx

% =========================================================================
% LOOP SETUP FOR THE THESIS
% =========================================================================
% We turn lista_p into a cell array so that we can mix numbers and
% text
lista_p = {1, 3, 5, 'Gauss'};
% 3 weights x 4 phi = 12 tables.
lista_weights = {'WendlandC2', 'WendlandC4', 'WendlandC6'};
nP = numel(lista_p);
nW = numel(lista_weights);

% One convergence table per combination: an nK x 6 matrix with
% [N, M, h, RMS, MAX, min(RCOND)].
Res = cell(nP, nW);%4x3 cell
for a = 1:numel(Res), Res{a} = zeros(nK, 6); end %we initialise each table
% with 0 in rows and columns.

warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix'); %we turn the warnings off.
fprintf('Iniciando simulaciones para el TFG...\n\n');

% 3. SWEEP IN N: the geometry that depends only on
% N is computed once per N (not 12 times). Inside each patch, the
% local interpolant is computed once per phi and reused for the
% weights.
for ik = 1:nK
    N = (2^k_list(ik) + 1)^s;    % Number of data points (Halton nodes)
n = max(round((N/N_j)^(1/s)), recinfo.n0);% centres per
% axis (total: n^s)
wep = n;
M = n^s;

% Data nodes using the Halton sequence
dsites = haltonseq(N, 2);
rhs = func_test_franke(dsites(:,1), dsites(:,2)); %evaluates Franke's
% function at the
% Halton nodes

% Centres of the patches
grid_pu = linspace(0, 1, n);
[xpu, ypu] = meshgrid(grid_pu);
cellctrs = [xpu(:), ypu(:)];
cellradius = 1/wep; % The radius of the patch is alpha*H=1/n.
rho = cellradius; % used in the normalisation (r/rho)^p.

% real fill distance h_{chi,Omega} (for the tables and the order)
h = distancia_llenado_knnsearch(dsites, epoints);

% --- 4. SPATIAL SEARCH ---
% Trees are created to find quickly which nodes fall inside
% each patch.
% The tree of nodes and DM_eval are built once
% per N (they used to be rebuilt for each combination).
datatree = createns(dsites, 'NSMethod', 'kdtree');
DM_eval = DistanceMatrixCSRBF(epoints, cellctrs, wep); % sparse
% matrix used to build the weights

% Only the patches that TOUCH some point of the grid matter; and for
% the data nodes of each patch we precompute their indices once.
patchcols = find(any(DM_eval, 1));%it checks column by column whether in each
% column (associated with a patch) there is some element different from zero,
% that is, some evaluation point in that patch. patchcols stores
% the indices (the columns) that satisfy the above.
patch_idx = cell(numel(patchcols), 1); %creates an array.
for jj = 1:numel(patchcols)
    j = patchcols(jj);%refers to the j-th centre
    [ic, ~] = rangesearch(datatree, cellctrs(j,:), cellradius);%gives us
    % the index of the Halton nodes that are at a distance smaller than
    % cellradius from the centre c_j.
    patch_idx{jj} = ic{1}';%puts into the array the indices of the Halton
    % points inside patch j.
    if isempty(patch_idx{jj}), DM_eval(:, j) = 0; end     % if the patch
    % has no Halton points, we set that column to zero, so that we do not
    % take it into account later in the MPU.
end

% --- 5. CONSTRUCTION OF THE WEIGHT MATRICES (SEM) ON epoints ---
% One SEM matrix per weight, computed once per N (it does not depend
% on phi). SEM{iw}(i,j) = omega_j(epoints(i,:)) (normalised Shepard
% weight).
SEM = cell(1, nW);%is an array with one row and as many columns as weights.
for iw = 1:nW
    weight_type = lista_weights{iw};

    % --- DEFINITION OF THE WEIGHT FUNCTIONS ---
    % The variable 'r' of these anonymous functions will directly receive
    % the sparse matrix returned by DistanceMatrixCSRBF.
    % That is, 'r' is u = 1 - r/rho.
    switch weight_type
        case 'WendlandC2'
            % Formula: r.^4 .* (5*spones(r) - 4*r)
            wf = @(ep, r) r.^4 .* (5*spones(r) - 4*r);

          case 'WendlandC4'
              % Adaptation for C4
              wf = @(ep, r) r.^6 .* (35*r.^2 - 88*r + 56*spones(r));

          case 'WendlandC6'
              % Adaptation for C6
              wf = @(ep, r) r.^8 .* (-32*r.^3 + 121*r.^2 - 154*r ...
                + 66*spones(r));
    end

    % It is evaluated by passing 'ep' and the matrix 'DM_eval'
    SEM_unscaled = wf(ep, DM_eval);
    S_inv = 1 ./ (SEM_unscaled * ones(M, 1));
    SEM{iw} = spdiags(S_inv, 0, neval^2, neval^2) * SEM_unscaled;
  %Matrix whose ij-th entry is (SEMiw)ij=omega_j(epoints(i,:)).
end

% We initialise the global interpolant vector (one per combination).
Pf = cell(nP, nW);%4x3 array, one Pf per weight and phi.
for a = 1:numel(Pf), Pf{a} = zeros(neval^2, 1); end
min_rcond = Inf(1, nP);% row vector of nP elements all equal to
% infinity.

% --- 6. MAIN LOOP OVER THE LOCAL PATCHES ---
%The local geometry (nodes, dist_j, evaluation points) is computed
% once per patch, and the local interpolant once per phi. Only the
% final assembly (with SEM{iw}) goes through the weights.
for jj = 1:numel(patchcols)
    j = patchcols(jj);%index of a patch with some evaluation point.

    % Find data inside cell j (precomputed)
    idx = patch_idx{jj}; % indices of Halton nodes in patch j.
    if isempty(idx), continue; end %if there are no nodes, skip what comes
    % below.

% Local data nodes and interpolation
local_dsites = dsites(idx, :);%coordinates of the nodes of patch
% j.
N_j_loc = length(idx);

% Matrix A_j of the local method
dist_j = DistanceMatrix(local_dsites, local_dsites);

% Find evaluation points falling in patch j (from the
% pattern of column j of DM_eval, which matches the support of
% SEM).
eidx = find(DM_eval(:, j));%gives the indices of the columns (patches) with
% some element different from zero. By what we did before, these
% are patches with some Halton node and some evaluation point.
if isempty(eidx), continue; end %if it is empty, go to the next patch
% with the for loop.
local_epoints = epoints(eidx, :);%evaluation points in the patch
% j.

% Matrix B (distances between evaluation points and data nodes)
dist_eval = DistanceMatrix(local_epoints, local_dsites);%sparse
% distance between nodes and evaluation points.

for ip = 1:nP
    p = lista_p{ip};%iterates over the list of weights.

    % Implementation of the Gaussian case and Phi=r^p.
    if isnumeric(p)%It tells apart the cases p=1,3,5 from the Gaussian one.
        % We now take Phi_rho=(r/rho)^p.
        A_j = (dist_j/rho).^p; % Phi(r) = (r/rho)^p
        % Matrix P_j defined for the new local basis of
        % polynomials.
        P_j = build_poly_basis(local_dsites, p, ...
          cellctrs(j,:), rho);
    else
        A_j = exp(-(ep * dist_j).^2); % Phi(r) = exp(-(ep*r)^2)
        % The Gaussian does not need a polynomial basis, empty matrix
        P_j = build_poly_basis(local_dsites, 0, cellctrs(j,:), ...
          rho);
    end

    L = size(P_j, 2); % The number of columns gives us the number
    % of basis vectors, hence the dimension of the space of
    % polynomials.

    % =========================================================
    % THEORETICAL CHECKS
    % ---------------------------------------------------------
    % To guarantee that the local system is non-singular, the
% set of nodes must be Pi_{m-1}^s-unisolvent.
% It only applies if there are polynomials (L > 0).
if L > 0
    % Condition 1: there must be at least as many points as
    % elements in the polynomial basis.
    if N_j_loc < L
        warning(['Parche %d: N_j (%d) < L (%d). ' ...
          'Sistema infradeterminado para los polinomios.'], ...
          j, N_j_loc, L);
    end
    % Condition 2: the matrix P_j must have full rank.
    % If rank(P_j) < L, it means that the local points lie
    % on an algebraic curve of degree m-1 and
    % are not unisolvent. The system will fail.
    if rank(P_j) < L
        warning(['Parche %d: Puntos locales NO ' ...
          'unisolventes (rank(P_j)=%d < L=%d).'], j, ...
          rank(P_j), L);
    end
end
% =========================================================

% Assembly of the local matrix
% IM = [A_j, P_j ; P_j^T, 0]
IM = [A_j, P_j; P_j', zeros(L, L)];

% Right-hand side (f_j extended with zeros for the
% polynomial constraints).
f_j = [rhs(idx); zeros(L, 1)];

% Conditioning diagnosis of the local system.
rc = rcond(full(IM));%full IM turns it from sparse into
% full.
if rc < min_rcond(ip), min_rcond(ip) = rc; end %for a fixed p,
%if rc is smaller than the previous one, we choose this rc as the worst.

% Solution of the local system: coefs=(lambda,c)'
coefs = IM \ f_j;

% --- Evaluation of the local interpolant at the evaluation
% points of the patch.
% Evaluation using the Gaussian Phi or Phi=r^p
if isnumeric(p)
    % With the same normalised Phi as in A_j (for consistency).
    B = (dist_eval/rho).^p; % Phi(r) = (r/rho)^p.
    % Matrix Q (polynomial basis evaluated at the evaluation
    % points)
    Q = build_poly_basis(local_epoints, p, cellctrs(j,:), rho);
else
    B = exp(-(ep * dist_eval).^2);
                        Q = build_poly_basis(local_epoints, 0, cellctrs(j,:), ...
                          rho); % empty Q
                  end

                  % Assembly of the evaluation matrix, that is, the
                  % matrix which, multiplied by coefs, gives us the
                  % local approximant j, evaluated at each of the
                  % evaluation points contained in patch j.
                  EM = [B, Q];

                  % Local interpolant
                  localfit = EM * coefs;

                  % Accumulation into the global interpolant using the
                  % Shepard weights (SEM). The same localfit works for the
                  % nW types of weights (it does not depend on the weight); only SEM{iw} changes.
                  for iw = 1:nW
                      Pf{ip,iw}(eidx) = Pf{ip,iw}(eidx) + localfit .* ...
                        SEM{iw}(eidx, j);%computes
                      % Global approximant on the evaluation points, taking
                      % into account the contributions of points falling in
                      % several patches.
                  end
            end
      end

      % --- 6. COMPUTATION OF ERRORS ---
      for ip = 1:nP
          for iw = 1:nW
              maxerr = norm(Pf{ip,iw} - exact, inf);
              rms_err = norm(Pf{ip,iw} - exact) / neval;%because
              % sqrt{neval^s)=neval,
              % for s=2
              Res{ip,iw}(ik, :) = [N, M, h, rms_err, maxerr, min_rcond(ip)];
          end
      end
      fprintf(' k=%d N=%d n=%d M=%d h=%.3e hecho\n', k_list(ik), N, ...
        n, M, h);
end

warning('on', 'MATLAB:nearlySingularMatrix');
warning('on', 'MATLAB:singularMatrix');

% =========================================================================
% PRINTING OF THE 12 TABLES. p_obs = log(RMS_prev/RMS)/log(h_prev/h).
% =========================================================================
for ip = 1:nP
    if isnumeric(lista_p{ip}), nombre_base = sprintf('r^%d', lista_p{ip});
    else, nombre_base = 'Gaussiana'; end
    for iw = 1:nW
            D = Res{ip, iw};
            pobs = NaN(nK,1);
            for r = 2:nK, pobs(r) = log2(D(r-1,4)/D(r,4)); end
            fprintf(['\n===============================================' ...
              '==================\n']);
            fprintf('TABLA: Base %s | Peso: %s\n', nombre_base, ...
              lista_weights{iw});
            fprintf(['================================================' ...
              '=================\n']);
            T = table(D(:,1), D(:,2), D(:,4), pobs, D(:,5), D(:,6), ...
                'VariableNames', {'N','M','RMS_error','orden_p', ...
                'max_error','min_RCOND'});
            disp(T);
      end
end

% =========================================================================
%LATEX OUTPUT FOR OVERLEAF.
% One table per phi, with the 3 weights (C2,C4,C6) in blocks of (RMS,p,max).
% p_obs = log(RMS_prev/RMS)/log(h_prev/h). Copy and paste directly.
% =========================================================================
fprintf('\n\n');
fprintf(['%%%%%%%%%%%%%%%% TABLAS EN LATEX (copiar a Overleaf) ' ...
  ' %%%%%%%%%%%%%%%%\n\n']);
pesos_tex = {'C^2', 'C^4', 'C^6'};
for ip = 1:nP
    if isnumeric(lista_p{ip})
        base_lbl = sprintf('r%d', lista_p{ip});
        cap_phi = sprintf('$\\Phi=r^{%d}$', lista_p{ip});
    else
        base_lbl = 'gauss';
        cap_phi = sprintf(['la gaussiana no estacionaria ($\\' ...
          'varepsilon=%.4f$)'], ep);
    end

      % Observed order per weight (one p column per block)
      Pobs = cell(1, nW);
      for iw = 1:nW
          D = Res{ip, iw};
          pv = NaN(nK, 1);
          for r = 2:nK, pv(r) = log2(D(r-1,4)/D(r,4)); end
          Pobs{iw} = pv;
      end

      fprintf('\\begin{table}[H]\n');
      fprintf(' \\centering\n');
      fprintf(' \\scalebox{0.68}{\n');
      fprintf(' \\begin{tabular}{| c | c | cccc | cccc | cccc |}\n');
      fprintf(' \\hline\n');
      fprintf([' & & \\multicolumn{4}{c|}{Wendland $%s$}' ...
               ' & \\multicolumn{4}{c|}{Wendland $%s$}' ...
               ' & \\multicolumn{4}{c|}{Wendland $%s$} \\\\\n'], ...
               pesos_tex{1}, pesos_tex{2}, pesos_tex{3});
      fprintf([' $k$ & $N$ & RMS-error & $p$ & max-error & RCOND' ...
               ' & RMS-error & $p$ & max-error & RCOND' ...
               ' & RMS-error & $p$ & max-error & RCOND \\\\\n']);
      fprintf(' \\hline\n');
      for ik = 1:nK
          fprintf(' %d & %d', k_list(ik), Res{ip,1}(ik,1));
          for iw = 1:nW
              D = Res{ip, iw};
              if ik == 1
                  fprintf(' & %.4e & & %.4e & %.2e', D(ik,4), D(ik,5), ...
                    D(ik,6));
              else
                  fprintf(' & %.4e & %.4f & %.4e & %.2e', D(ik,4), ...
                    Pobs{iw}(ik), D(ik,5), D(ik,6));
              end
          end
          fprintf(' \\\\\n');
      end
      fprintf(' \\hline\n');
      fprintf(' \\end{tabular}}\n');
      fprintf([' \\caption{RMSE y error m\\''aximo frente al ' ...
        'n\\''umero de nodos de ' ...
               'Halton $N$ para el MPU con %s y pesos de Wendland $C^2$,' ...
               ' $C^4$ y $C^6$, ' ...
               'sobre la funci\\''on de Franke.}\n'], cap_phi);
      fprintf(' \\label{tab:mpu_franke_%s}\n', base_lbl);
      fprintf('\\end{table}\n\n');
end

% =========================================================================
% GLOBAL RBF (comparison): Gaussian with fixed eps and r^p, p=1,3,5.
% Following code 10 of the appendix: direct global RBF interpolation, without
% polynomial augmentation (for r^p we solve DM.^p \ rhs, like the distance
% matrix with r^1). Its own sweep up to an affordable N_max (dense system).
% The eps of the Gaussian is fixed as the smallest one with rcond >= rcond_target at
% N_max (rcond grows with eps and N_max is the worst case: this way rcond > threshold over
% the whole sweep). p_obs = log2(RMS_prev/RMS).
% =========================================================================
k_rbf     = 1:7;          % N=(2^k+1)^2 => N=9..16641 (affordable N_max)
nK_rbf    = numel(k_rbf);
lista_rbf = {'Gauss', 1, 3, 5};
nP_rbf    = numel(lista_rbf);

% FIXED eps of the Gaussian: smallest eps with rcond >= rcond_target at N_max.
N_mx_rbf = (2^max(k_rbf) + 1)^s;
nd_mx_rbf = haltonseq(N_mx_rbf, s);
DM_mx_rbf = DistanceMatrix(nd_mx_rbf, nd_mx_rbf);
grc = @(e) rcond(exp(-(e*DM_mx_rbf).^2));
e_lo = ep_native; e_hi = 2*ep_native;
while grc(e_hi) < rcond_target, e_hi = 2*e_hi; if e_hi > 1e4, break;
end, end
for it = 1:50
    em = 0.5*(e_lo + e_hi);
    if grc(em) >= rcond_target, e_hi = em; else, e_lo = em; end
end
eps_fijo = e_hi;
fprintf(['\neps_fijo (RBF global gaussiano) = %.4f (rcond en' ...
  ' N_max=%d: %.2e)\n'], ...
  eps_fijo, N_mx_rbf, grc(eps_fijo));

% Computation and storage: Res_rbf{ip} = [k, N, RMS, MAX, RCOND]
Res_rbf = cell(1, nP_rbf);
for ip = 1:nP_rbf
    phi_rbf = lista_rbf{ip};
    Drbf = zeros(nK_rbf, 5);
    for ik = 1:nK_rbf
        N = (2^k_rbf(ik) + 1)^s;
        nodes = haltonseq(N, s);
        rhs_rbf = func_test_franke(nodes(:,1), nodes(:,2));
        DM_rbf = DistanceMatrix(nodes, nodes);
        DM_ev   = DistanceMatrix(epoints, nodes);
        if isnumeric(phi_rbf)
            IM_rbf = DM_rbf.^phi_rbf;
            Pf_rbf = (DM_ev.^phi_rbf) * (IM_rbf \ rhs_rbf);
        else
            IM_rbf = exp(-(eps_fijo*DM_rbf).^2);
            Pf_rbf = exp(-(eps_fijo*DM_ev).^2) * (IM_rbf \ rhs_rbf);
        end
        rms_rbf = norm(Pf_rbf - exact) / neval;
        max_rbf = norm(Pf_rbf - exact, inf);
        Drbf(ik, :) = [k_rbf(ik), N, rms_rbf, max_rbf, rcond(IM_rbf)];
    end
    Res_rbf{ip} = Drbf;
end

% Text tables (console)
for ip = 1:nP_rbf
    if isnumeric(lista_rbf{ip}), nombre_rbf = sprintf('r^%d', ...
        lista_rbf{ip});
    else, nombre_rbf = 'Gaussiana'; end
    D = Res_rbf{ip};
    pobs = NaN(nK_rbf, 1);
    for rr = 2:nK_rbf, pobs(rr) = log2(D(rr-1,3)/D(rr,3)); end
    fprintf(['\n===============================================' ...
      '==================\n']);
    fprintf('TABLA: RBF GLOBAL | Base %s\n', nombre_rbf);
    fprintf(['================================================' ...
        '=================\n']);
      T = table(D(:,1), D(:,2), D(:,3), pobs, D(:,4), D(:,5), ...
          'VariableNames', {'k','N','RMS_error','orden_p','max_error' ...
          ,'RCOND'});
      disp(T);
end

% LaTeX tables (to copy into Overleaf), one per Phi, with an RCOND column
fprintf('\n\n');
fprintf(['%%%%%%%%%%%%%%%% TABLAS RBF GLOBAL EN LATEX (copiar' ...
  ' a Overleaf) %%%%%%%%%%%%%%%%\n\n']);
for ip = 1:nP_rbf
    if isnumeric(lista_rbf{ip})
        base_lbl = sprintf('r%d', lista_rbf{ip});
        cap_phi = sprintf('$\\Phi=r^{%d}$', lista_rbf{ip});
    else
        base_lbl = 'gauss';
        cap_phi = sprintf('la gaussiana ($\\varepsilon=%.4f$)', eps_fijo);
    end
    D = Res_rbf{ip};
    pv = NaN(nK_rbf, 1);
    for r = 2:nK_rbf, pv(r) = log2(D(r-1,3)/D(r,3)); end
    fprintf('\\begin{table}[H]\n');
    fprintf(' \\centering\n');
    fprintf(' \\scalebox{0.9}{\n');
    fprintf(' \\begin{tabular}{| c | c | c c c c |}\n');
    fprintf(' \\hline\n');
    fprintf(' $k$ & $N$ & RMS-error & $p$ & max-error & RCOND \\\\\n');
    fprintf(' \\hline\n');
    for ik = 1:nK_rbf
        if ik == 1
            fprintf(' %d & %d & %.4e & & %.4e & %.2e \\\\\n', ...
              D(ik,1), D(ik,2), D(ik,3), D(ik,4), D(ik,5));
        else
            fprintf(' %d & %d & %.4e & %.4f & %.4e & %.2e \\\\\n', ...
              D(ik,1), D(ik,2), D(ik,3), pv(ik), D(ik,4), D(ik,5));
        end
    end
    fprintf(' \\hline\n');
    fprintf(' \\end{tabular}}\n');
    fprintf([' \\caption{RMSE y error m\\''aximo frente al n' ...
      '\\''umero de nodos de ' ...
             'Halton $N$ para el RBF global con %s, sobre la func' ...
             'i\\''on de Franke.}\n'], cap_phi);
    fprintf(' \\label{tab:rbf_franke_%s}\n', base_lbl);
    fprintf('\\end{table}\n\n');
end

% =========================================================================
% TIMING STUDY (tic/toc): MPU (Wendland C6) vs global RBF, per Phi.
% One full solve per (N, Phi) is timed [geometry + systems +
% evaluation], which is what sets the real cost. The median of
% several repetitions is taken. The number of repetitions is not fixed per N (the
% times of MPU and RBF are not comparable with each other): each case is repeated
% until a budget t_budget is used up or max_rep is reached, so that a
%fast solve (MPU) is averaged over many repetitions and a slow one (RBF)
%over few (often 1), automatically. The first sample of each case
%is a warm-up and is discarded.
%Note: within each method the times of the 4 Phi come out almost equal (the
% MPU is dominated by the geometry, common to all 4; the RBF by the dense NxN solve,
%also common). The 4 columns serve as a check of that independence from Phi.
%Besides, there are more MPU rows than RBF ones: the MPU sweeps all of k_list while
% the RBF is only timed up to N_rbf_cap (affordable dense system).
% =========================================================================
k_tim      = k_list;                  % same sweep in k as the MPU
nK_tim     = numel(k_tim);
Nvec_tim = (2.^k_tim + 1).^s;         % N per row
N_rbf_cap = 16641;                    % RBF only up to an affordable N (dense NxN)
t_budget = 0.5;                      % time budget per case [s]
max_rep    = 15;               % cap on repetitions (small and fast N)
Tmpu = NaN(nK_tim, nP);            % MPU C6 times [s], one column per Phi
Trbf = NaN(nK_tim, nP);        % global RBF times [s], one column per Phi
fprintf('\nCronometrando MPU (C6) y RBF global...\n');
for ik = 1:nK_tim
    N = Nvec_tim(ik);
    for ip = 1:nP
         p = lista_p{ip};
         % ---------------- MPU (Wendland C6), Phi = p ----------------
         smp = zeros(1, max_rep); ns = 0; tacc = 0; warm = true;
         while true %the loop repeats indefinitely until a break that comes
           % further down.
             tt = tic;
             n = max(round((N/N_j)^(1/s)), recinfo.n0); M = n^s; rho = 1/n;
             dsites = haltonseq(N, s);
             rhs = func_test_franke(dsites(:,1), dsites(:,2));
             grid_pu = linspace(0, 1, n); [xpu, ypu] = meshgrid(grid_pu);
             cellctrs = [xpu(:), ypu(:)];
             datatree = createns(dsites, 'NSMethod', 'kdtree');
             DM_eval = DistanceMatrixCSRBF(epoints, cellctrs, n);
             SEM_u = DM_eval.^8 .* (-32*DM_eval.^3 + 121*DM_eval.^2 ...
                - 154*DM_eval + 66*spones(DM_eval));
             S_inv = 1 ./ (SEM_u * ones(M, 1));
             SEM = spdiags(S_inv, 0, neval^2, neval^2) * SEM_u;
             Pf = zeros(neval^2, 1);
             patchcols = find(any(DM_eval, 1));
             for jj = 1:numel(patchcols)
                  j = patchcols(jj);
                  ic = rangesearch(datatree, cellctrs(j,:), rho);
                  idx = ic{1}'; if isempty(idx), continue; end
                  local_dsites = dsites(idx, :);
        dist_j = DistanceMatrix(local_dsites, local_dsites);
        eidx = find(DM_eval(:, j)); if isempty(eidx), continue; end
        local_epoints = epoints(eidx, :);
        dist_eval = DistanceMatrix(local_epoints, local_dsites);
        if isnumeric(p)
             A_j = (dist_j/rho).^p;
             P_j = build_poly_basis(local_dsites, p, ...
               cellctrs(j,:), rho);
             B = (dist_eval/rho).^p;
             Q = build_poly_basis(local_epoints, p, ...
               cellctrs(j,:), rho);
        else
             A_j = exp(-(ep*dist_j).^2);
             P_j = build_poly_basis(local_dsites, 0, ...
               cellctrs(j,:), rho);
             B = exp(-(ep*dist_eval).^2);
             Q = build_poly_basis(local_epoints, 0, ...
               cellctrs(j,:), rho);
        end
        L = size(P_j, 2);
        coefs = [A_j, P_j; P_j', zeros(L, L)] \ [rhs(idx);
          zeros(L, 1)];
        localfit = [B, Q] * coefs;
        Pf(eidx) = Pf(eidx) + localfit .* SEM(eidx, j);
    end
    e = toc(tt);
    if warm
        warm = false;             % discards the 1st one (warm-up)
    else
        ns = ns + 1; smp(ns) = e; tacc = tacc + e;
        if tacc >= t_budget || ns >= max_rep, break; end
    end
end
Tmpu(ik, ip) = median(smp(1:ns));
% ---------------- global RBF, Phi = p ----------------
if N <= N_rbf_cap
    smp = zeros(1, max_rep); ns = 0; tacc = 0; warm = true;
    while true
        tt = tic;
        nodes = haltonseq(N, s);
        rhs_r = func_test_franke(nodes(:,1), nodes(:,2));
        DMd = DistanceMatrix(nodes, nodes);
        DMe = DistanceMatrix(epoints, nodes);
        if isnumeric(p)
             Pfr = (DMe.^p) * (DMd.^p \ rhs_r);
        else
             Pfr = exp(-(eps_fijo*DMe).^2) * (exp(-(eps_fijo* ...
               DMd).^2) \ rhs_r);
        end
        e = toc(tt);
                     if warm
                          warm = false;
                     else
                          ns = ns + 1; smp(ns) = e; tacc = tacc + e;
                          if tacc >= t_budget || ns >= max_rep, break; end
                     end
                  end
                  Trbf(ik, ip) = median(smp(1:ns));
            end
      end
end

% --- normalised ratios: MPU by N*log(N); RBF by N^3 ---
Rmpu = NaN(nK_tim, nP); Rrbf = NaN(nK_tim, nP);
for ik = 2:nK_tim
    fN = (Nvec_tim(ik)*log(Nvec_tim(ik))) / (Nvec_tim(ik-1)*log ...
      (Nvec_tim(ik-1)));
    fN3 = Nvec_tim(ik)^3 / Nvec_tim(ik-1)^3;
    for ip = 1:nP
        Rmpu(ik,ip) = (Tmpu(ik,ip)/Tmpu(ik-1,ip)) / fN;
        Rrbf(ik,ip) = (Trbf(ik,ip)/Trbf(ik-1,ip)) / fN3;
    end
end

% --- TABLE 1 (LaTeX): times [s] ---
fprintf('\n\n');
fprintf(['%%%%%%%%%%%%%%%% TABLAS DE TIEMPOS EN LATEX (copiar a' ...
  ' Overleaf) %%%%%%%%%%%%%%%%\n\n']);
fprintf('\\begin{table}[H]\n \\centering\n \\scalebox{0.8}{\n');
fprintf(' \\begin{tabular}{| c | c | cccc | cccc |}\n \\hline\n');
fprintf([' & & \\multicolumn{4}{c|}{MPU ($C^6$) [s]}' ...
         ' & \\multicolumn{4}{c|}{RBF global [s]} \\\\\n']);
fprintf([' $k$ & $N$ & $r^1$ & $r^3$ & $r^5$ & Gauss' ...
         ' & $r^1$ & $r^3$ & $r^5$ & Gauss \\\\\n']);
fprintf(' \\hline\n');
for ik = 1:nK_tim
    rowstr = sprintf(' %d & %d', k_tim(ik), Nvec_tim(ik));
    for ip = 1:nP, rowstr = [rowstr, sprintf(' & %.3e', Tmpu(ik,ip))]; end
    for ip = 1:nP
        if isnan(Trbf(ik,ip)), rowstr = [rowstr, ' & '];
        else, rowstr = [rowstr, sprintf(' & %.3e', Trbf(ik,ip))]; end
    end
    fprintf('%s \\\\\n', rowstr);
end
fprintf(' \\hline\n \\end{tabular}}\n');
fprintf([' \\caption{Tiempos de computo (tic/toc) del MPU (Wendland' ...
  ' $C^6$) y' ...
         ' del RBF global, por $\\Phi$, sobre la funcion de Franke.}\n']);
fprintf(' \\label{tab:tiempos}\n\\end{table}\n\n');

% --- TABLE 2 (LaTeX): normalised ratios (they should tend to ~1) ---
fprintf('\\begin{table}[H]\n \\centering\n \\scalebox{0.8}{\n');
fprintf(' \\begin{tabular}{| c | c | cccc | cccc |}\n \\hline\n');
fprintf([' & & \\multicolumn{4}{c|}{MPU: $(t_{k}/t_{k-1})/' ...
  '(N\\log N)$-ratio}' ...
         ' & \\multicolumn{4}{c|}{RBF: $(t_{k}/t_{k-1})/N^3$-ratio}' ...
         ' \\\\\n']);
fprintf([' $k$ & $N$ & $r^1$ & $r^3$ & $r^5$ & Gauss' ...
         ' & $r^1$ & $r^3$ & $r^5$ & Gauss \\\\\n']);
fprintf(' \\hline\n');
for ik = 1:nK_tim
    rowstr = sprintf(' %d & %d', k_tim(ik), Nvec_tim(ik));
    for ip = 1:nP
        if isnan(Rmpu(ik,ip)), rowstr = [rowstr, ' & '];
        else, rowstr = [rowstr, sprintf(' & %.4f', Rmpu(ik,ip))]; end
    end
    for ip = 1:nP
        if isnan(Rrbf(ik,ip)), rowstr = [rowstr, ' & '];
        else, rowstr = [rowstr, sprintf(' & %.4f', Rrbf(ik,ip))]; end
    end
    fprintf('%s \\\\\n', rowstr);
end
fprintf(' \\hline\n \\end{tabular}}\n');
fprintf([' \\caption{Cocientes de tiempos normalizados: MPU por' ...
  ' $N\\log N$ y' ...
         ' RBF global por $N^3$. Valores $\\approx 1$ confirman el' ...
         ' orden esperado.}\n']);
fprintf(' \\label{tab:cocientes}\n\\end{table}\n\n');
