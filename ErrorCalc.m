% This function computes the errors based on all permutations of the true answer. 
% A bipartite matching helps to compute misclassified nodes.

function [bip_score]=ErrorCalc(vecA, vecB, m)
   cluster_mat = zeros(m);

for i=1:m    
    set_true = find(vecA==i);
    for j=1:m
        set_ans = find(vecB == j);
        common_indx = intersect(set_ans,set_true);
        cluster_mat(i,j) = length(common_indx);
    end
end

[val,~,~]=bipartite_matching(cluster_mat);
bip_score=val/sum(cluster_mat(:));


end