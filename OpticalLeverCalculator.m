global L w0x w0y lamda omegaM Gamma meff rho a b h bw P WtoV km k T kb hbar c SthFF Responsitivity PDGain
L = 200e-3;
w0x = 60e-6;
w0y = 60e-6;
lamda = 1064e-9;
omegaM = 182194*2*pi;
Gamma = 5;
rho = 3100;
a = 1.5e-3;
b = 3.5e-3;
h = 100e-9;
bw =1;
P = 1730e-6;
Responsitivity = 0.5;
PDGain = 300e3;
WtoV = Responsitivity*PDGain;
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
ds = 0:wx/20:3*wx;
peakHeightC = zeros(1, length(ds));
noiFloorC = zeros(1, length(ds));
samplepoints = [omegaM - 2*pi*5*bw:2*pi*bw:omegaM + 2*pi*5*bw];

%{
samplepoints = [omegaM/(2*pi) - 1500*bw:bw:omegaM/(2*pi) + 1500*bw];
d = 0;
SPF = F(0,0,1,0,d);
GeomNoise = GeomNoi(d, SPF);
PSD_theory = arrayfun(@(x) S(x,pi/2, SPF, GeomNoise), 2*pi*samplepoints);
%}

for i = 1 : length(ds)
    d = ds(i);
    SPF = F(0,0,1,0,d);
    GeomNoise = GeomNoi(d, SPF);
    heights_C = arrayfun(@(x) S(x,pi/2, SPF, GeomNoise), samplepoints);
    power = sum(10.^(heights_C/10))/20;
    peakHeightC(i) = power;
    noiFloorC(i) = 10^(S(omegaM - 100*Gamma, pi/2, SPF, GeomNoise)/10)/20;
    toc 
end
plot(noiFloorC*0.14./peakHeightC);
ds = ds/w(w0x, L);
clear ans d f GeomNoise re SNR SPF wx y DN
clear global
%plot(freq(check1:check2-1), PSD(check1:check2-1), samplepoints-8.515, PSD_theory')

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
    global P bw k km w0x SthFF WtoV hbar c Responsitivity
    re = 10*log10(20) + ...
    10*log10(2*bw*((2*P*SPF*k*km*w0x*abs(chi(omega)).*sin(theta)).^2*SthFF/2 +...
    P*SPF^2*hbar*c*k*(abs(2*P/(c*k)*(k*km*w0x)^2*chi(omega).*sin(theta) - cos(theta)).^2/Responsitivity + sin(theta).^2/Responsitivity + GeomNoise/Responsitivity))*(WtoV)^2); %in unit of dB/bin
end

function re = Sth(f, theta, SPF, GeomNoise)
    global P bw k km w0x SthFF WtoV hbar c
    re = 10*log10(20*bw*4*P^2*SPF^2*(k*km*w0x)^2*SthFF*(abs(chi(2*pi*f))).^2*(WtoV)^2);
end