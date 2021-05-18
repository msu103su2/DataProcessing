function [AC, DC] = CM_Cali(Fweight, dy, w0y, bw, Pac, Pdc)
    A = 2*SPF*dy/w0y*cos(theta);
    AC = 0.5*(Pac*A)^2/bw;
    DC = (Pdc*A)^2;
end

function re = u(m,n,x,y)
    wx = w(w0x, L);
    wy = w(w0y, L);
    re = (2/pi)^(1/2)/sqrt(2^(m+n)*factorial(m)*factorial(n)*wx*wy)*exp(-x.^2/wx^2-y.^2/wy^2).*hermiteH(m,sqrt(2)*x/wx).*hermiteH(n,sqrt(2)*y/wy);
end

function re = F(m,n,l,k,d, doff)
    global w0x w0y L
    wx = w(w0x, L);
    wy = w(w0y, L);
    scale = 3;
    func = @(x, y) u(m,l,x,x).*Fweight(x, d, doff);
    re = integral(func, -scale*wx, scale*wx,'RelTol',1e-3);
end

function re = Fweight(x, d, doff)
    re = heaviside(x-doff - d/2) - heaviside(-(x-doff + d/2));
end