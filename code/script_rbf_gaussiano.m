% RBFInterpolation2D

clearvars; close all; clc;

%% 1. SETUP AND DEFINITIONS
% Parameters of the RBF (Gaussian)
rbf_gaussiana = @(e,r) exp(-(e*r).^2);
eps_fijo = 21.1; % Shape parameter (epsilon)

% Franke test function
f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
func_test_franke = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);

resol_malla = 40;
valores_k_graficas = 1:6;
valores_k_tabla1 = 1:7;
resol_epsilon_vec = 250; % High resolution in order to capture peaks of
% instability

% Unified 3D camera views
view_surfaces_3d = [-40, 30];
view_points_3d = [-40, 30];

fprintf('\n=== INICIANDO HERRAMIENTA DE INTERPOLACIÓN RBF 2D ===\n');

%% ================= PART 1: GRID GENERATION ========
epoints = MakeSDGrid(2, resol_malla);
% We prepare the matrices xe, ye for the plotting functions.
% IMPORTANT: since MakeSDGrid gives us in each row the coordinates (x,y)
% of the points of the equispaced grid, just to make the plots
% MATLAB needs the data in matrix form in order to know in which order
% they go
xe = reshape(epoints(:,1), resol_malla, resol_malla);
ye = reshape(epoints(:,2), resol_malla, resol_malla);
% We evaluate the exact solution (test function) on the fine evaluation
% grid
exact = func_test_franke(epoints(:,1), epoints(:,2));


%% ================= PART 2: 3D SURFACES ========
fprintf(['\n--- PARTE 2: Generando Gráficas de Superficie y Puntos de' ...
  ' Halton ---\n']);
% Loop for k from 1 to 6
for i = 1:length(valores_k_graficas)
    k = valores_k_graficas(i);
    N = (2^k+1)^2;
    % Interpolation nodes
    dsites = haltonseq(N, 2);
    ctrs = dsites;

    if k == 1 || k == 2
        figure('Name', sprintf('Nodos y Plano z=0 (k=%d)', k), ...
          'NumberTitle', 'off');
        % 1. We draw the transparent plane at z=0
        surf(xe, ye, zeros(size(xe)), 'FaceAlpha', 0.2, 'EdgeColor', ...
          'none', 'FaceColor', 'interp');
        hold on; %So as to draw on top
         % 2. We draw the points in 3D fixing z=0
        plot3(dsites(:,1), dsites(:,2), zeros(N, 1), 'ko', ...
          'MarkerFaceColor', 'b', 'MarkerSize', 6);
        axis([0 1 0 1 -0.2 1.2]);
        grid on;
            view(view_points_3d);
            xlabel('x'); ylabel('y'); zlabel('z');
            title(sprintf('Nodos de Halton y Plano z=0 (k = %d, N = %d)', ...
              k, N));
            set(gca, 'FontSize', 12);
            hold off;
      end

      rhs = func_test_franke(dsites(:,1), dsites(:,2)); %test function
      % evaluated at the Halton points
      % We compute the distance matrices only once for both methods
      DM_data = DistanceMatrix(dsites, ctrs);
      DM_eval = DistanceMatrix(epoints, ctrs);

      % --- METHOD 1: Gaussian ---
      IM_gauss = rbf_gaussiana(eps_fijo, DM_data); %matrix of the method
      % (evaluated at the nodes) where IM_ij is the Gaussian evaluated at
      % ||x_i-x_j||
      EM_gauss = rbf_gaussiana(eps_fijo, DM_eval); %matrix where EM_ij is the
      % Gaussian evaluated at ||z_i-x_j||, with z_i the points of the equispaced
      % grid
      warning('off', 'all');
      Pf_gauss = EM_gauss * (IM_gauss \ rhs); %gives us the value of the
      % Gaussian RBF interpolating function evaluated at the points of the equispaced
      % grid, to compare them later with the values of the original
      % test function (Franke's function)
      warning('on', 'all');
      maxerr_gauss = norm(Pf_gauss - exact, inf);
      caption_gauss = sprintf('Gaussiana (\\epsilon = 21.1) N = %d', N);
      PlotSurf(xe, ye, Pf_gauss, exact, maxerr_gauss, view_surfaces_3d, ...
        caption_gauss);

      % --- METHOD 2: Distance matrix ---
      Pf_dm = DM_eval * (DM_data \ rhs); %gives us the value of the
      % interpolating function of distance matrices evaluated at the points of the
      % equispaced grid, to compare them later with the values of the
      % original test function (Franke's function)
      maxerr_dm = norm(Pf_dm - exact, inf);
      caption_dm = sprintf('Matriz Distancia con N = %d', N);
      PlotSurf(xe, ye, Pf_dm, exact, maxerr_dm, view_surfaces_3d, ...
        caption_dm);
end


%% ================= PART 3: FIRST LATEX TABLE ========
fprintf(['\n\n--- PARTE 3: Generando Tabla LaTeX 2.1 (Calculando k=7, ' ...
  'tardará un poco) ---\n']);
warning('off', 'all');

fprintf('\n%% === TABLA 2.1 ===\n\n');
fprintf('\\begin{table}[H]\n');
fprintf(' \\centering\n');
fprintf(' \\scalebox{0.9}{\n');
fprintf(' \\begin{tabular}{cc | ccc | ccc}\n');
fprintf(' \\hline\n');
fprintf([' & & \\multicolumn{3}{c|}{Gaussiano} & \\multicolumn{3}{c}' ...
  '{Matrices distancia} \\\\\n']);
fprintf([' $k$ & $N$ & RMS-error & $p$ & max-error & RMS-error &' ...
  ' $p$ & max-error \\\\\n']);
fprintf(' \\hline\n');

M = resol_malla^2;

for k_tab = valores_k_tabla1
    N_tab = (2^k_tab + 1)^2;
    dsites_tab = haltonseq(N_tab, 2);
    ctrs_tab = dsites_tab;
    %Same process as before, but now going up to k=7, in order to make the table
    rhs_tab = func_test_franke(dsites_tab(:,1), dsites_tab(:,2));

   DM_data_tab = DistanceMatrix(dsites_tab, ctrs_tab);
   DM_eval_tab = DistanceMatrix(epoints, ctrs_tab);

   IM_gauss = rbf_gaussiana(eps_fijo, DM_data_tab);
   EM_gauss = rbf_gaussiana(eps_fijo, DM_eval_tab);
   Pf_gauss = EM_gauss * (IM_gauss \ rhs_tab);

   rms_gauss = norm(Pf_gauss - exact) / sqrt(M);
   max_gauss = norm(Pf_gauss - exact, inf);

   Pf_dm = DM_eval_tab * (DM_data_tab \ rhs_tab);
   rms_dm = norm(Pf_dm - exact) / sqrt(M);
   max_dm = norm(Pf_dm - exact, inf);

   % COMPUTATION AND PRINTING WITH CONVERGENCE RATE 'p'
   if k_tab == valores_k_tabla1(1)
        % First row: without 'p' (blank gap)
        fprintf(' %d & %d & %1.6e & & %1.6e & %1.6e & & %1.6e \\\\\n', ...
                 k_tab, N_tab, rms_gauss, max_gauss, rms_dm, max_dm);
   else
        % Computation of 'p' using the previous RMS
        p_gauss = log(rms_gauss / rms_gauss_prev) / log(0.5);
        p_dm = log(rms_dm / rms_dm_prev) / log(0.5);
        fprintf([' %d & %d & %1.6e & %1.4f & %1.6e & %1.6e & %1.4f & ' ...
          '%1.6e \\\\\n'], ...
                 k_tab, N_tab, rms_gauss, p_gauss, max_gauss, rms_dm, ...
                 p_dm, max_dm);
   end

   % We store the values for the next iteration (for the first one they
    % are not needed, because there is no element in the p column
    rms_gauss_prev = rms_gauss;
    rms_dm_prev = rms_dm;
end
fprintf(' \\hline\n \\end{tabular}}\n\\end{table}\n');
warning('on', 'all');


%% ================= PART 4: ERROR VS EPSILON CURVES ========
fprintf('\n--- PARTE 4: Generando Curvas Error vs Epsilon ---\n');
fprintf('    (Evaluando matrices, por favor espera un momento...)\n');

eps_vec = linspace(0.1, 25, resol_epsilon_vec); %I do not start at 0, because
% there IM_ij is 1 for all its elements, so it is non-singular.
%I have split the interval into 250 points to have plenty of precision when
%observing peaks of instability in the error-epsilon plots and when
%making the table

figure('Name', 'Errores vs Epsilon (Escalado Dinamico)', 'Position', ...
  [100, 100, 1200, 700]);

for i = 1:length(valores_k_graficas)
    k = valores_k_graficas(i);
    N = (2^k+1)^2;
    fprintf('Calculando curva para N = %d ... ', N);

    dsites = haltonseq(N, 2);
    rhs = func_test_franke(dsites(:,1), dsites(:,2));
    DM_data = DistanceMatrix(dsites, dsites);
    DM_eval = DistanceMatrix(epoints, dsites);

    rms_vec = zeros(size(eps_vec));
    max_vec = zeros(size(eps_vec));

    warning('off', 'all');
    for j = 1:length(eps_vec)
        ep_val = eps_vec(j);
        IM = rbf_gaussiana(ep_val, DM_data);%Matrix in which each component
        % is the Gaussian with the epsilon the for loop is at at that
        % moment, evaluated at the elements of the distance matrix of the method
        % (evaluated at the Halton points)
        EM = rbf_gaussiana(ep_val, DM_eval); %Matrix in which each component
        % is the Gaussian with the epsilon the for loop is at at that
        % moment, evaluated at the elements of the distance matrix evaluated at the points
        % equispaced grid
        Pf = EM * (IM \ rhs); %Interpolating function evaluated at the points
        % equispaced grid

        rms_vec(j) = norm(Pf - exact) / sqrt(M);
        max_vec(j) = norm(Pf - exact, inf);
      end
      fprintf('Hecho.\n');
      warning('on', 'all');

      subplot(2, 3, i);%plot with 2 rows and 3 columns of graphs. i indicates
      % the position of each graph as the loop runs from 1 to 6
      semilogy(eps_vec, max_vec, 'k--', 'LineWidth', 1.5); hold on;
      semilogy(eps_vec, rms_vec, 'k-', 'LineWidth', 1.5);
      %we have superimposed the maximum error and the RMS error on one plot, where
      %the y axis is in logarithmic scale
      title(sprintf('N = %d', N));
      xlabel('\epsilon'); ylabel('error');


      y_min = min(rms_vec);
      y_max = max(max_vec);
      %Adjust the limit of the plots
      if y_min <= 0 || isnan(y_min) || isinf(y_min), y_min = 1e-16; end
      if isnan(y_max) || isinf(y_max), y_max = 1e2; end

      xlim([0 25]);
      ylim([y_min / 5, y_max * 5]); % Fully dynamic adjustment guaranteed
      % -----------------------------------

      if i == 1
          legend('max-error', 'RMS-error', 'Location', 'southwest'); %makes
          % the legend of the solid and dashed lines appear only in the first
          % plot (remove if I want)
      end
      grid on;
end


%% ================= PART 5: SECOND LATEX TABLE ========
fprintf('\n\n--- PARTE 5: Generando Tabla LaTeX 2.2 ---\n');

fprintf('\n%% === TABLA 2.2 ===\n\n');
fprintf('\\begin{table}[H]\n');
fprintf(' \\centering\n');
fprintf(' \\scalebox{0.85}{\n');
fprintf(' \\begin{tabular}{cc | cccc | cccc}\n');
fprintf(' \\hline\n');
fprintf([' & & \\multicolumn{4}{c|}{$\\epsilon$ ``seguro'' más' ...
  ' pequeño} & \\multicolumn{4}{c}{smallest RMS-error $\\epsilon$}' ...
  ' \\\\\n']);
fprintf([' $k$ & $N$ & $\\epsilon$ & RMS-error & $p$ & max-error &' ...
  ' $\\epsilon$ & RMS-error & $p$ & max-error \\\\\n']);
fprintf(' \\hline\n');

for i = 1:length(valores_k_graficas)
    k_tab = valores_k_graficas(i);
    N_tab = (2^k_tab + 1)^2;
    dsites_tab = haltonseq(N_tab, 2);
    ctrs_tab = dsites_tab;
    rhs_tab = func_test_franke(dsites_tab(:,1), dsites_tab(:,2));
    DM_data_tab = DistanceMatrix(dsites_tab, ctrs_tab);
    DM_eval_tab = DistanceMatrix(epoints, ctrs_tab);

    warning('off', 'all');
    %-----------Looking for the minimum epsilon
    options = optimset('Display','off', 'TolX', 1e-3); %Display off is
    % so that MATLAB does not report every step it takes. TolX 10^(-3) is the
    % precision: when you are so sure of where the minimum is that the doubt
    % is smaller than 0.001, it stops searching.
    rms_obj = @(e) norm((rbf_gaussiana(e, DM_eval_tab) * ...
      (rbf_gaussiana(e, DM_data_tab) \ rhs_tab)) - exact) / sqrt(M);%Everything
    %that we had done before in several steps it does in a single
    %formula. We give it the value of epsilon and it computes the RMS error
    eps_rms_min = fminbnd(rms_obj, 0.05, 30, options);%It looks for the minimum of
    %the function of the previous line between 0.05 and 30 with options, which is what
    %we had told it before about the tolerance and the Display off

    IM_opt = rbf_gaussiana(eps_rms_min, DM_data_tab);
    EM_opt = rbf_gaussiana(eps_rms_min, DM_eval_tab);
    Pf_opt = EM_opt * (IM_opt \ rhs_tab);
    rms_opt = norm(Pf_opt - exact) / sqrt(M);%we compute the RMS error for
    % the minimum epsilon
    max_opt = norm(Pf_opt - exact, inf);
    warning('on', 'all');

    eps_safe = 25;
    for test_ep = 25:-0.01:0.01%we start looking for epsilon at 25
      % (very stable) and we go down
         IM_test = rbf_gaussiana(test_ep, DM_data_tab);
         if rcond(IM_test) < 1e-13 %if rcond is smaller than 10^(-13) (we have
           % chosen this one to be safe, even though the precision is 10^(-16)),
           % we keep the epsilon just above it
             eps_safe = test_ep + 0.01;
             break;
         end
    end
%if we get to the end of the loop and rcond has never been smaller than the
% fixed value, eps_safe is still 25, but really the safest one
% should be 0.01; we take 0.02 to be safe, and we see that in no case
% do we reach that limit
    if eps_safe == 25 && rcond(rbf_gaussiana(0.01, DM_data_tab)) > 1e-13
         eps_safe = 0.02;
    end

    warning('off', 'all');
      IM_safe = rbf_gaussiana(eps_safe, DM_data_tab);
      EM_safe = rbf_gaussiana(eps_safe, DM_eval_tab);
      Pf_safe = EM_safe * (IM_safe \ rhs_tab);
      %we compute the errors for the safe epsilon
      rms_safe = norm(Pf_safe - exact) / sqrt(M);
      max_safe = norm(Pf_safe - exact, inf);
      warning('on', 'all');

      % COMPUTATION AND PRINTING WITH CONVERGENCE RATE 'p'
      if k_tab == valores_k_graficas(1)
           % Row 1: without 'p' (blank gap)
           fprintf([' %d & %d & %.2f & %1.6e & & %1.6e & %.2f & %1.6e' ...
             ' & & %1.6e \\\\\n'], ...
                    k_tab, N_tab, eps_safe, rms_safe, max_safe, ...
                    eps_rms_min, rms_opt, max_opt);
      else
           % Computation of 'p' using the previous RMS
           p_safe = log(rms_safe / rms_safe_prev) / log(0.5);
           p_opt = log(rms_opt / rms_opt_prev) / log(0.5);
           fprintf([' %d & %d & %.2f & %1.6e & %1.4f & %1.6e & %.2f &' ...
             ' %1.6e & %1.4f & %1.6e \\\\\n'], ...
                    k_tab, N_tab, eps_safe, rms_safe, p_safe, max_safe, ...
                    eps_rms_min, rms_opt, p_opt, max_opt);
      end

      % We store the values for the next iteration
      rms_safe_prev = rms_safe;
      rms_opt_prev = rms_opt;
end

fprintf(' \\hline\n');
fprintf(' \\end{tabular}}\n');
fprintf(['\\caption{Interpolación por RBF Gaussiano a la función de' ...
  ' Frankee.}\n']);
fprintf('\\end{table}\n');
