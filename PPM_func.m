%% Main module of the projected power method for joint alignment proposed in the paper
% ``The Projected Power Method: An Efficient Algorithm for Joint Alignment from Pairwise Differences'' by Y. Chen and E. J. Candes.
% The initialization method is modified so that the algorithm by Candes can
% accommodate community detection in censored block model.

function zt = PPM_func(m, n, param, X_gt, L_opt, proj, option, X0_in)
%  Projected power method for joint alignment
%    param:   parameters
%    X_gt:    ground truth (matrix representation)
%    L_opt:   multiplication operator ( L_opt(x) = L* x )
%    proj:    projection operator 

%%  Spectral initialization
% Compute rank-m approximation of L

if  (strcmp(option,'QR')==1)
    U = randn(n*m, m);
    B = L_opt( U );

    % Block power method with the assistance of QR decomposition
    for iter = 1: param.max_poweriter
        [U,Sigma]  = qr(B, 0);     % economy-size QR decomposition
        B = L_opt( U );

        if norm(B - U*Sigma) / norm(Sigma) <= param.eps_power, break;  end
    end
   X0 = U * diag(diag(Sigma)) * U( randi(n*m), :)'; 
   
elseif(strcmp(option,'Rand') ==1)
    
    I_m = eye(m);
    for i = 1:n    
        X0((i-1)*m+1:i*m,:) = I_m(:,randi([1,m]));
    end

elseif(strcmp(option,'EigenMax') ==1)
    
     I_mn = eye(m*n);
     B = L_opt(I_mn);
     [X0,~] = eigs(B,1);

elseif(strcmp(option,'Spectral') ==1)
    
    I_m = eye(m);
    X0 = zeros(n*m,1);
    for i =1:n
       X0((i-1)*m+1:i*m) = I_m(:,X0_in(i));
    end 
    
end
% Round a random column to generate an initial guess

zt = proj( X0 ); 

%% Projected power method
for iter = 2:  (param.maxIter - 1)
    zold = zt;      
    zt   = L_opt( zt );     
    zt   = proj( zt );  
    if norm(zold - zt) == 0,  break;  end
end

% Round the final iterate to ensure feasibility

zt = ProjectOperator( zt, n, m, param.mu0 );


