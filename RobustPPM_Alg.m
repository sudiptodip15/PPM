function RobustPPM_Alg(Problem, param)

    L  = Problem.L;
    Y = Problem.Y;   
    m = param.m;
    n = param.n;
    outlier_index = param.outlier_index;
        
    Y_true = Y;
    L_true = L;  
    true_node = 1:n;
    adv_node = [];    
    z_ans = zeros(param.n,1);
    
    %% Iterative Refinement       
    max_iter=param.max_robust_iter;
    for iter=1:max_iter
        
        % Spectral initialization
        
        [Y_tilde,~] = eigs(Y_true,m); 
        [IDX,~] = kmeans(Y_tilde,m);

         % Define PPM function handles

        Lambda      = svds( L_true, 2);      % Compute 2nd largest singularvalue of L
        L_opt       = @(I) L_true  * I;      % Matrix-vector multiplication operator
        proj        = @(z) ProjectOperator(z, n, m, param.mu0/Lambda(2)); % Projection operator
        
        % Run projected power method
        init = 'Spectral';
        X_init = IDX;
        z  = PPM_func(m, n, param, Problem.X_gt, L_opt, proj, init, X_init);  % z is the output and relErrors is the L2 error

        
        for i = 1:n
            temp = z((i-1)*m+1:i*m);
            [~,z_ans(true_node(i))]=max(temp);
        end
        
        if(iter == 1)    
            bip_score=ErrorCalc(Problem.x_gt(1:outlier_index), z_ans(1:outlier_index), m);
            fprintf('PPM without Robustness , Error = %f ', (1-bip_score));
        end
        
        if(iter ~= 1)
            comp_mat_true = repmat(z_ans(true_node),1,param.m);
            comp_mat_cluster = repmat(1:param.m,n,1);
            index_same = comp_mat_true == comp_mat_cluster;
            index_diff = comp_mat_true ~= comp_mat_cluster;
            for i = 1 : length(adv_node)            
                 mask = repmat(Y(true_node,adv_node(i)),1,param.m);
                 score = sum(mask.*index_same - mask.*index_diff);
                 [~,z_ans(adv_node(i))] = max(score);             
            end
            
             bip_score=ErrorCalc(Problem.x_gt(1:param.outlier_index), z_ans(1:param.outlier_index), m);
             fprintf('Iteration %d , Error = %f ', iter-1, (1-bip_score));
        
        end
               
       
        % Compute Likelihood Score

        comp_mat_1 = repmat(z_ans',param.n,1);
        comp_mat_2 = repmat(z_ans,1,param.n);
        index_same = comp_mat_1 == comp_mat_2;   
        index_diff = comp_mat_1 ~= comp_mat_2;

        mask_same = Y.*index_same;     
        mask_diff = Y.*index_diff;
        if(strcmp(param.adv_model,'rnd')==1)
            like_score = sum(mask_same) - sum(mask_diff);
        else            
            like_score = sum(mask_same==1);
        end
        
        [~,sort_index]=sort(like_score,'descend');
        true_node = sort_index(1:outlier_index);
        adv_node = sort_index(outlier_index+1:param.n);   
        
        false_positive= sum(adv_node<=outlier_index);
        fprintf('  False positive = (%d / %d )\n', false_positive, (param.n-outlier_index));
        
        %disp(adv_node);
        
        % Retain rows and columns of true_nodes from Y_true and L_true 
        
        Y_true = Y(true_node,true_node);  
        n=length(true_node);
        L_true = zeros(n*m,n*m);
        for i=1:n
            for j=i+1:n
                L_true((i-1)*m+1:i*m,(j-1)*m+1:j*m) = L((true_node(i)-1)*m+1:true_node(i)*m,(true_node(j)-1)*m+1:true_node(j)*m);
            end
        end
        
        L_true = sparse(L_true);
        L_true = L_true+L_true';       
        
        
    end   
    fprintf(' Iterative Rebustness Algorithm Reduces Error !!!\n');

end