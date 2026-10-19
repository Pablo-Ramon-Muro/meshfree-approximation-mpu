% PlotError2D(xe,ye,Pf,exact,maxerr,fview,caption)
% Produces a plot of the absolute error for the surface Pf
% (abs(Pf-exact))
function PlotError2D(xe, ye, Pf, exact, maxerr, fview, caption)

      figure;

      % Dynamic reshape adapted to the grid xe ---
      Error_2D = reshape(abs(Pf - exact), size(xe));

      errorplot = surf(xe, ye, Error_2D);
      set(errorplot, 'FaceColor', 'interp', 'EdgeColor', 'none');

      % Compatibility: we keep caxis ---
      [cmin, cmax] = caxis;
      caxis([cmin - 0.25 * maxerr, cmax]);

      view(fview);
      colormap hsv;

      vcb = colorbar('vert');

      % Avoids an error if maxerr is 0 ---
      if maxerr > 0
          ylim(vcb, [0, maxerr]);
      end

      set(get(vcb, 'YLabel'), 'String', 'Error', 'FontSize', 14);

      title(caption);
end
