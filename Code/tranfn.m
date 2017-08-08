% Create second characteristic function approximation
% Returns H matrix (*1.43)

function [N, D] = tranfn(A, p)
    syms s N D
 
    N    = cell(p+1,1);
    N{1} = A{2,1};
    N{2} = A{2,1}*A{2,1};
    
    D    = cell(p+1,1);
    D{1} = A{1,1};
    D{2} = A{1,1}*A{2,1} + s*A{3,1};
    
    for i = 3:p
        N{i} = A{i,1}*N{i-1} + s*A{i+1,1}*N{i-2};
        D{i} = A{i,1}*D{i-1} + s*A{i+1,1}*D{i-2};
    end
end
