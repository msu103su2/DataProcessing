function re = fit_logLorentzian(x, y, fmI)
    f = @(fm, gamma, offset, x) logLorentzian(x, fm, gamma, offset);
    re = fit(x, y, fittype(f), 'StartPoint', [fmI, 20, 180], 'Lower', [fmI-50, 0, 120], 'Upper', [fmI+50,300, 280]);
end

function re = logLorentzian(x, fm, gamma, offset)
    re = offset - 10*log10((fm^2 - x.^2).^2*(2*pi)^4+(2*pi*x*gamma).^2);
end