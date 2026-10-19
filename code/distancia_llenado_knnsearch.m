function h = distancia_llenado_knnsearch(chi, Y)
% h_{chi,Omega} = max_i min_j ||Y(i,:) - chi(j,:)|| computed very
% efficiently
%
% chi: N x s, nodes
% Y:   Neval x s, evaluation points (fine grid approximating Omega)
% Requires: Statistics and Machine Learning Toolbox.

      [~, D] = knnsearch(chi, Y);     % D(i) = distance from Y(i,:) to its
      % nearest node in chi
      h = max(D);
end
