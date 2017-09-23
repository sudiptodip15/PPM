function w = ProjectOperator(z, n, m, mu)
%  blockwise projection onto standard simplex
%    z:  a vector to be projected
%    mu: scaling factor

w = zeros(n*m,1);

% project onto standard simplex
for irows = 1: n
    rows = ((irows-1)* m + 1) : (irows * m) ;
    w(rows) = ProjectOntoSimplex( mu* z(rows) );
end 

