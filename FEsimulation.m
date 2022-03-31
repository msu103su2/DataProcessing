function unew = FEsimulation(HG, z, t, n, h, ht, x, y)
    k = 2*pi/1.064e-6;
    ui = HG.FEu(x, y, z, k);
    uold = ui;
    unew = ui;
    figure;
    set(gcf,'position',[500 0 1200 1000]);
    delta =  2*h*del2(uold, ht)/(1i*n*k);
    %s = surf(abs(unew));
    %s.XData = x;
    %s.YData = y;
    surf(abs(unew));
    hold on
    while z < t
        delta = 2*h*del2(uold, ht)/(1i*n*k);
        unew = uold + delta;
        z = z + h;
        disp([max(abs(uold), [], 'all'), max(abs(delta), [], 'all'), z]);
        uold = unew;
        surf(abs(unew));
        drawnow();
        pause(0.01);
    end
    if t == 0
        unew = ui;
    end
end