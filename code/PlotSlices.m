function PlotSlices(xe,ye,ze,Pf,xslice,yslice,zslice,caption)
% PlotSlices(xe,ye,ze,Pf,neval,xslice,yslice,zslice,caption)
% Produces a volume slice plot for Pf
% xslice,yslice,zslice define the range and the number of slices.

figure
% Using size(xe) makes the reshape failure-proof.
Pf_3D = reshape(Pf, size(xe));

pfit = slice(xe, ye, ze, Pf_3D, xslice, yslice, zslice);
set(pfit, 'FaceColor', 'interp', 'EdgeColor', 'none')
daspect([1 1 1])
view(3);
axis([0 1 0 1 0 1])

vcb = colorbar;
ylabel(vcb, 'Function value', 'FontSize', 14);
title(caption)

end
