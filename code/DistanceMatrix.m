function DM = DistanceMatrix(dsites,ctrs)
[M,s] = size(dsites); [N,s] = size(ctrs);
DM = zeros(M,N);
% Builds the distance matrix, where the dsites are the points where the
%function ||x-x_k|| is evaluated, while the centres are the x_k. Both usually
%coincide
% The ndgrid command produces two matrices of size MxN:
%   dr: made of N identical columns (each one holds
%        the d-th coordinate of the M data points)
%   cc: made of M identical rows (each one holds
%        the d-th coordinate of the N centres)
for d=1:s
    [dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
    DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);
