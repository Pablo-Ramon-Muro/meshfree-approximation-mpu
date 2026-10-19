function gridpoints = MakeSDGrid(s, neval)
% MakeSDGrid(s,neval)
% Produces a matrix of equispaced points in the cube of dimension s (one
% point per row)

    if (s == 1)
        gridpoints = linspace(0, 1, neval)';
        return;
    end

    % 1. We create a 1D base vector
    v = linspace(0, 1, neval);

    % 2. We preallocate a 'Cell Array' to collect the 's' outputs of
    % ndgrid
    % This replaces creating [x1, x2, ..., xs]
    grid_cells = cell(1, s);

    % 3. We call ndgrid expanding the vector 'v' s times as input,
      % and we collect the output in the cells.
      % The {:} syntax expands the cell automatically.
      inputs = repmat({v}, 1, s);
      [grid_cells{:}] = ndgrid(inputs{:});

      % 4. We preallocate the final matrix of points for maximum speed
      gridpoints = zeros(neval^s, s);

      % 5. We dump the data by columns into the final matrix
      for d = 1:s
          gridpoints(:, d) = grid_cells{d}(:);
      end

end
