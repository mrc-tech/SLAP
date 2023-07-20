clear
clc


conjgrad([4,1;1,3], [1;2], [2;1]) % funziona
% conjgrad([1,3,3;4,5,6;7,8,9], [1;1;1], [0;0;0]) % NON FUNZIONA




function x = conjgrad(A, b, x)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;

        if(mod(i,50) == 0); r = b - A*x; else
            r = r - alpha * Ap; end
        rsnew = r' * r;
        
        if sqrt(rsnew) < 1e-10
            break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end