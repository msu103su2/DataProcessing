function re = fit_logLorentzian(x, y, fmI)
    f = @(fm, gamma, offset, x) logLorentzian(x, fm, gamma, offset);
    re = fit(x, y, fittype(f), 'StartPoint', [fmI, 20, 180], 'Lower', [fmI-50, 0, 0], 'Upper', [fmI+50,300, 280]);
    f = @(fm, gamma, A, noi, x) logLorentzian_noi(x, fm, gamma, A, noi);
    re = fit(x, y, fittype(f), 'StartPoint', [fmI, 1, 1e10, 5e-9], 'Lower', [fmI-50, 0, 1, 1e-25], 'Upper', [fmI+50,300, 1e28, 1]);
end

function re = logLorentzian(x, fm, gamma, offset)
    re = offset - 10*log10((fm^2 - x.^2).^2*(2*pi)^4+(2*pi*x*gamma).^2);
end

function re = logLorentzian_noi(x, fm, gamma, A, noi)
    re =  10*log10(A./((fm^2 - x.^2).^2*(2*pi)^4+(2*pi*x*gamma).^2)+noi);
end