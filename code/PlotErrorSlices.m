% PlotErrorSlices(xe,ye,ze,Pf,exact,xslice,yslice,zslice,caption)
% Produces a volumetric error slice plot: abs(Pf-exact)
function PlotErrorSlices(xe, ye, ze, Pf, exact, xslice, yslice, ...
  zslice, caption)

      % Slice plot for the error
      figure;

      % We compute the absolute error
      error_abs = abs(Pf - exact);

      % We use size(xe) to avoid a dimension error
      error_3D = reshape(error_abs, size(xe));

      % We create the slices
      errorplot = slice(xe, ye, ze, error_3D, xslice, yslice, zslice);
      set(errorplot, 'FaceColor', 'interp', 'EdgeColor', 'none');

      daspect([1 1 1]);
      view(3);
      axis([0 1 0 1 0 1]);

      % --- We go back to the caxis syntax ---
      [cmin, cmax] = caxis;
      caxis([cmin - 0.25 * cmax, cmax]);

      colormap hsv;

      % Colour bar
      vcb = colorbar('vert');

      % We only apply ylim if cmax is greater than 0, to avoid internal
      % failures
      if cmax > 0
          ylim(vcb, [0, cmax]);
      end

      set(get(vcb, 'YLabel'), 'String', 'Error', 'FontSize', 14);
      title(caption);

end
