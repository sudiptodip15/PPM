function [Problem] = ProbGenAdverse(param)
%  Generate synthetic examples
%    param:  parameters
%    Problem: returns the synthetic example

n = param.n;        
m = param.m;
outlier_index = param.outlier_index;
Problem.param   = param;

% Generate the ground truth
Problem.x_gt    = randi( m, [n,1] ); 
Problem.x_gt(outlier_index+1:n)=m+1;

% Matrix/vector representation of the truth
for i = 1: outlier_index 
    rowIdx = ((i-1)*m + 1) : (i*m);
    Problem.X_gt( rowIdx, : ) = circshift( eye(m),  Problem.x_gt(i)); 
end

for i = outlier_index+1: n
    rowIdx = ((i-1)*m + 1) : (i*m);
    Problem.X_gt( rowIdx, : ) = ones(m);
end

% Compute the input data matrix L
[Problem.L, Problem.Y] = generate_input(Problem.x_gt, param);




function [L0,Y] = generate_input(x_gt, param)
% Generate the input matrix L (as a sparse matrix)
%   x_gt:  ground truth
%   param: input parameters

n       = param.n;  
m       = param.m;  
outlier_index = param.outlier_index;
p_obs   = param.p_obs;
P0      = param.P0 / sum( param.P0 );
Pcs     = cumsum(P0);   

   
% Omega is the index set of the sampled pairs
Omega = zeros(n,n);
for i = 1: n
    for j = (i+1): n
        if rand() < p_obs,  Omega(i,j) = 1; end
    end
end

[rowsOmega, colsOmega] = find(Omega);
nSamples   = length(rowsOmega);
    
multi_label_PPM = 0;

% When multi_label_PPM = 1, the original PPM formulation is implemented.
% Setting multi_label_PPM = 0 allows the censored block model to work.

if(multi_label_PPM == 1)
    L0entries  = ones( nSamples*m, 2);
    
    % Store the indices of nonzero entries for all cyclic permutation matrices
    Eye_all_idx = zeros(m,2,m);
    for i = 1: m
       [row_idx, col_idx] = find( circshift( eye(m), i-1) ); 
       Eye_all_idx(:, 1, i) = row_idx;
       Eye_all_idx(:, 2, i) = col_idx;
    end

    % Generate all non-zero entries of L
    for l = 1: nSamples
        row_idx = ((l-1)*m+1) : (l*m);
        off  = find( Pcs > rand() , 1 ) - 1 ; 
        meas = mod( x_gt(rowsOmega(l)) - x_gt(colsOmega(l)) + off, m);
        
        % the starting points of the block (minus 1 in each dimensino)
        row_min = ( rowsOmega(l) -1 ) * m; 
        col_min = ( colsOmega(l) -1 ) * m;    
        L0entries( row_idx, 1:2 ) = repmat([row_min, col_min], [m,1]) + Eye_all_idx(:, :, meas + 1);
    end

    % Generate L0 as a sparse matrix
    L0 = sparse(L0entries(:,1), L0entries(:,2), ones(nSamples*m,1), n*m, n*m);
    
    % Matrix Y is un-information in this case
    Y = eye(n,n);
    fprintf('Y matrix in uninformative. Spectral methods not to be used !\n');
    
else
      
    debias = 1;       % '1' subtracts off the mean from likelihood matrix.
    L0 = zeros(n*m,n*m);    
    Y = zeros(n,n);    
    
    like_pos = eye(m);
    like_neg = ones(m)-eye(m);    

    if(debias == 1)
       like_pos = like_pos - (1/m^2)*(sum(like_pos(:)))*ones(m);
       like_neg = like_neg - (1/m^2)*(sum(like_neg(:)))*ones(m);    
    end
    
    Eye_two_idx(:,:,1) = like_pos;
    Eye_two_idx(:,:,2) = like_neg;
    
    % Generate all non-zero entries of L
    
    
    for l = 1: nSamples        
        off  = find( Pcs > rand() , 1 ) - 1 ;    
        
       % Only two options, either offset is 0, or non-zero
       if(off~=0); off = 1; end;                   

       if(rowsOmega(l) > outlier_index || colsOmega(l) > outlier_index)
           if(strcmp(param.adv_model,'rnd')==1)
                meas=randi([0,1]);
           else
                meas=mod(off+1,2);
           end
       else
           meas = mod((x_gt(rowsOmega(l))~= x_gt(colsOmega(l))) + off,2);  
       end
       
       Y(rowsOmega(l),colsOmega(l))= 1 - 2*meas;
       L0(((rowsOmega(l)-1)* m+1):(rowsOmega(l)*m), ((colsOmega(l)-1) * m+1):(colsOmega(l)* m)) = Eye_two_idx(:, :, meas + 1);
             
       
    end   
    
    Y=Y+Y'+eye(n);      
    L0 = sparse(L0);
end

L0 = L0 + L0';
