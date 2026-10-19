% Example script for chapter 1, with a for loop for the plots and for generating
% the LaTeX table
% -------------------------------------------------------------------------
clearvars; close all; clc;

%% ================= PART 1: GENERATION OF THE PLOTS =================
% We define the combinations of s, k and neval that we want to evaluate
valores_s = [1, 2, 3];
valores_k = [2, 4, 3];
valores_neval = [1000, 40, 10];

% Loop going through the 3 combinations in order to plot
for i = 1:length(valores_s)

      s = valores_s(i);
      k = valores_k(i);
      neval = valores_neval(i);

      fprintf(['\n[GRÁFICAS] Evaluando caso: s = %d, k = %d, neval' ...
        ' = %d\n'], s, k, neval);

    N = (2^k+1)^s;
    M = neval^s;

    dsites = haltonseq(N,s);%the nodes for the interpolation (Halton
    % points)
    ctrs = dsites;
    epoints = MakeSDGrid(s,neval);%equispaced grid in the cube
    rhs = testfunction(s,dsites);%exact values of the test function at the nodes

    IM = DistanceMatrix(dsites,ctrs);%matrix of the method, where IM_ij is
    % ||x_i-x_j||
    EM = DistanceMatrix(epoints,ctrs);%matrix where EM_ij is
% ||z_i-x_j||, with z_i the points of the equispaced grid


    Pf = EM * (IM\rhs); %gives us the value of the interpolating function evaluated at the
% points of the equispaced grid, to compare them later with the values
% of the original test function
    exact = testfunction(s,epoints);%test function evaluated at the points of the
    % equispaced grid

    maxerr = norm(Pf-exact,inf);
    rms_err = norm(Pf-exact)/sqrt(M);
    fprintf('RMS error:       %e\n', rms_err)
    fprintf('Error máximo: %e\n', maxerr)

    switch s
        case 1
            figure;
            plot(epoints, Pf)
            title(['Interpolación (s=', num2str(s), ', k=', ...
              num2str(k), ')']);
            figure;
            plot(epoints, abs(Pf-exact))
            title(['Error absoluto (s=', num2str(s), ', k=', ...
              num2str(k), ')']);
        case 2
            fview = [-30,30];
            xe = reshape(epoints(:,2),neval,neval);
            ye = reshape(epoints(:,1),neval,neval);
            caption = ['Interpolación (s=', num2str(s), ') false' ...
              ' colored by error.'];
            PlotSurf(xe, ye, Pf, exact, maxerr, fview, caption);
            caption = ['Error absoluto en la interpolación' ...
              ' (s=', num2str(s), ').'];
            PlotError2D(xe,ye,Pf,exact,maxerr,fview,caption);
        case 3
            xe = reshape(epoints(:,2),neval,neval,neval);
            ye = reshape(epoints(:,1),neval,neval,neval);
                ze = reshape(epoints(:,3),neval,neval,neval);
                xslice = .25:.25:1; yslice = 1; zslice = [0,0.5];
                caption = ['Slice plot de la interpolación' ...
                  ' (s=', num2str(s), ').'];
                PlotSlices(xe,ye,ze,Pf,xslice,yslice,zslice,caption);
                caption = ['Slice plot del error absoluto en la' ...
                  ' interpolación(s=', num2str(s), ').'];
                PlotErrorSlices(xe,ye,ze,Pf,exact,xslice,yslice, ...
                  zslice,caption);
            otherwise
                disp('-> No se muestran gráficos para s>3')
      end
end


%% ================= PART 2: GENERATION OF THE LATEX TABLE ============


% We stop MATLAB from filling the console with warnings about badly
% conditioned matrices
warning('off', 'all');

% Maximum limits of 'k' so as not to saturate the memory
k_max_por_dimension = [12, 6, 3, 2, 2, 1];

% Adjusted 'neval': fixed at 4 for s=4, 5 and 6
neval_tabla = [1000, 40, 10, 4, 4, 4];

fprintf('\\begin{table}[H]\n');
fprintf(' \\centering\n');
fprintf(' \\scalebox{0.8}{\n');
fprintf(' \\begin{tabular}{| c | cc | cc | cc | cc | cc | cc |}\n');
fprintf(' \\hline\n');
fprintf([' & \\multicolumn{2}{c|}{1D} & \\multicolumn{2}{c|}{2D}' ...
  ' & \\multicolumn{2}{c|}{3D} & \\multicolumn{2}{c|}{4D} & ' ...
  '\\multicolumn{2}{c|}{5D} & \\multicolumn{2}{c|}{6D} \\\\\n']);
fprintf([' $k$ & $N$ & RMS-error & $N$ & RMS-error & $N$ & RMS-error' ...
  ' & $N$ & RMS-error & $N$ & RMS-error & $N$ & RMS-error \\\\\n']);
fprintf(' \\hline\n');

% Loop for the rows (k from 1 to 12)
for k_tab = 1:12
    fprintf(' %2d ', k_tab); % Prints the value of k

      % Loop for the columns (dimensions from 1 to 6)
      for s_tab = 1:6
          if k_tab <= k_max_por_dimension(s_tab)
              N_tab = (2^k_tab + 1)^s_tab;
              neval_t = neval_tabla(s_tab);
              M_tab = neval_t^s_tab;

              dsites_tab = haltonseq(N_tab, s_tab);
              epoints_tab = MakeSDGrid(s_tab, neval_t);
              rhs_tab = testfunction(s_tab, dsites_tab);

              IM_tab = DistanceMatrix(dsites_tab, dsites_tab);
              EM_tab = DistanceMatrix(epoints_tab, dsites_tab);

              Pf_tab = EM_tab * (IM_tab \ rhs_tab);
              exact_tab = testfunction(s_tab, epoints_tab);

              rms_err_tab = norm(Pf_tab - exact_tab) / sqrt(M_tab);

              % Prints the values separating N and Error with the ampersand (&)
              fprintf('& %d & %1.6e ', N_tab, rms_err_tab);
          else
              % If we go past the limit, we leave the LaTeX cells empty
              % while respecting the columns
              fprintf('& & ');
          end
      end
      fprintf(' \\\\\\\\\n'); % Prints the LaTeX line break (\\)
end

fprintf(' \\hline\n');
fprintf(' \\end{tabular}}\n');
fprintf(['\\caption{Fit de la matriz distancia para los puntos de' ...
  ' Halton para dimensiones de 1 a 6.}\n']);
fprintf('\\end{table}\n');

warning('on', 'all');
