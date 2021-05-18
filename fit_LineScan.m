function re = fit_LineScan(step, power, m, n, a, b)
    f = @(stepSize, x0, k, s, A, x) fitFunc(x, stepSize, x0, k, s, A, m, n, a, b);
    re = fit(step, power, fittype(f), 'StartPoint', [0.8e-7, -5e-6, 0, b/2, 57], 'Lower', [50e-9, -a, -0.1, 0, 0], 'Upper', [1e-3, a, 0.1, b, 240]);
    plot(step, re(step), step, power, '*');
end

function power = fitFunc(step, stepSize, x0, k, s, A, m, n, a, b)
    theta = atan(k);
    x = step*stepSize - x0;
    y = k*x+s;
    power = cos(theta)*m*pi/a*cos(m*pi/a*x).*sin(n*pi/b*y)+sin(theta)*n*pi/b*sin(m*pi/a*x).*cos(n*pi/b*y);
    power = 10*log10(power.^2)+2*A;
end