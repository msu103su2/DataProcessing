function re = erfFit(x, power, Aangle, tiltangle, Pmax)
    x = x*1e-3;
    f = @(wx,P,x0,x) func(x, Aangle, tiltangle, wx, P, x0);
    re = fit(x, power', fittype(f), 'StartPoint', [1.25e-3, Pmax, max(x)], 'Lower', [3e-4, 0.9*Pmax, (max(x)-3)], 'Upper', [5e-3, 1.1*Pmax, (max(x)+3)]);
    fitx = 1e-3:1e-4:15e-3;
    plot(fitx,re(fitx),x,power,'*');
end

function re = func(x, Aangle, tiltangle, wx, P, x0)
    re = 1/2*P*(1+erf((x-x0)*sin(Aangle)*cos(tiltangle)/(wx/sqrt(2))));
end