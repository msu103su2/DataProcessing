global L w0x w0y lamda omegaM Gamma meff rho a b h bw P WtoV km k T kb hbar c SthFF
L = 55e-3;
w0x = 15e-6;
w0y = 15e-6;
lamda = 1064e-9;
omegaM = 122037*2*pi;
Gamma = 135.45;
rho = 3100;
a = 1.5e-3;
b = 3.5e-3;
h = 100e-9;
bw = 8.515;
P = 2037e-6;
WtoV = 42.7e3;
T = 300;
kb = 1.38064852e-23;
hbar = 1.054571817e-34;
c = 299792458;

meff = rho*a*b*h/4;
k = 2*pi/lamda;
km = 2*pi/a;
SthFF = 4*kb*T*Gamma*meff;

tic
wx = w(w0x,L);
xc = 0:wx/20:3*wx;
SNR = [];
peakHeightC = [];
noiFloorC = [];
DN = -110;
samplepoints = [omegaM - 5*bw:bw:omegaM + 5*bw];
for d = 0:wx/20:3*wx
    SPF = F(0,0,1,0,d);
    GeomNoise = GeomNoi(d, SPF);
    heights = arrayfun(@(x) S(x,pi/2, SPF, GeomNoise), samplepoints);
    power = 10*log10(sum(10.^(heights/10)));
    peakHeightC = [peakHeightC, power];
    noiFloorC = [noiFloorC, S(omegaM - 100*Gamma, pi/2, SPF, GeomNoise)];
    toc 
end
%peakHeightC = 10*log10(10.^(peakHeightC/10)+10^(DN/10));
%noiFloorC = 10*log10(10.^(noiFloorC/10)+10^(DN/10));
clear ans d f GeomNoise re SNR SPF wx y DN
clear global

function re = SNRc(omega, theta, d)
    global Gamma
    re = S(oemga, theta, d) - S(oemga-100*Gamma, theta, d);
end

function re = chi(omega)
    global meff omegaM Gamma
    re = 1/meff*((omegaM^2 - omega.^2)+1i*Gamma*omega).^(-1);
end

function re = GeomNoi(d, denominator)
    re = 0;
    for i = 1:5
        re = re+(abs(F(0,0,2*i+1,0,d)/denominator)).^2;
    end
end

function re = F(m,n,l,k,d)
    global w0x w0y L
    wx = w(w0x, L);
    wy = w(w0y, L);
    scale = 3;
    func = @(x, y) u(m,l,x,x).*Fweight(x,d);
    re = integral(func, -scale*wx, scale*wx,'RelTol',1e-3);
end

function re = Fweight(x,d)
    re = heaviside(x - d/2) - heaviside(-(x + d/2));
end

function re = u(m,n,x,y)
    global w0x w0y L
    wx = w(w0x, L);
    wy = w(w0y, L);
    re = (2/pi)^(1/2)/sqrt(2^(m+n)*factorial(m)*factorial(n)*wx*wy)*exp(-x.^2/wx^2-y.^2/wy^2).*hermiteH(m,sqrt(2)*x/wx).*hermiteH(n,sqrt(2)*y/wy);
end

function re = w(w0, z)
    global lamda
    re = w0*sqrt(1+(z*lamda/(pi*w0^2))^2);
end

function re = S(omega, theta, SPF, GeomNoise)
    global P bw k km w0x SthFF WtoV hbar c
    re = 10*log10(20) + ...
    10*log10(2*bw*((2*P*SPF*k*km*w0x*abs(chi(omega)).*sin(theta)).^2*SthFF +...
    P*SPF^2*hbar*c*k*(abs(2*P/(c*k)*(k*km*w0x)^2*chi(omega).*sin(theta) - cos(theta)).^2 + sin(theta).^2 + GeomNoise))*(WtoV)^2);
end