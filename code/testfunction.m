function tf = testfunction(s, points)
% tf = testfunction(s,points)
% Evaluates the test function (normalised, so that its maximum is 1) of the
% distance matrix example at the given points (points) for the
% given dimension (s)

tf = 4^s * prod(points .* (1 - points), 2);

end
