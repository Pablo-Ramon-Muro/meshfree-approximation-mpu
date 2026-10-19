function PlotSurf(xe, ye, Pf, exact, maxerr, fview, caption)
    figure
    % We use size(xe) so that it is robust
    Pf_plot = reshape(Pf, size(xe));
    Err_plot = reshape(abs(Pf - exact), size(xe));

      h = surf(xe, ye, Pf_plot, Err_plot);
      set(h, 'FaceColor', 'interp', 'EdgeColor', 'none');
      view(fview); colormap hsv;

      [cmin, cmax] = caxis;
      caxis([cmin - 0.25*maxerr, cmax]);
      vcb = colorbar;
      if maxerr > 0, ylim(vcb, [0, maxerr]); end
      ylabel(vcb, 'Error abs');
      title(caption);
end
