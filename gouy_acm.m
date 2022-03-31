function theta = gouy_acm(z, w0, lambda, f, d)
    T1 = [[1, lambda*d];[0, 1]];
    L1 = [[1, 0];[-1/(lambda*f), 1]];
    T2 = [[1, lambda*(z-d)];[0, 1]];
    z0 = pi*w0^2/lambda;
    if z < d
        theta = atan(z/z0);
    elseif z == d
        T = L1*T1;
        A = T(1,1);
        B = T(1,2);
        theta = atan(B/(A*pi*w0^2));
    else
        T = T2*L1*T1;
        A = T(1,1);
        B = T(1,2);
        theta = atan(B/(A*pi*w0^2));
    end
    if theta < 0
        theta = theta + pi;
    end
end