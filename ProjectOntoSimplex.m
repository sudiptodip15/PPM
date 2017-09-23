function w = ProjectOntoSimplex(v)
% Euclidean projection of a point v onto the standard simplex. 
% This is adapted from the code implemented by John Duchi

% Make sure all entries of v are positive by enforcing a global offset
v       = v - min(v) + 1;
u       = sort(v,'descend');
sv      = cumsum(u);
rho     = find(u > (sv - 1) ./ (1:length(u))', 1, 'last');
theta   = max(0, (sv(rho) - 1) / rho);
w       = max(v - theta, 0);
