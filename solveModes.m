function s = solveModes(fs)
    syms c a b;
    eqns = [fs(1) == c/(2*pi)*sqrt((pi/a)^2+(pi/b)^2),...
        fs(2) == c/(2*pi)*sqrt((2*pi/a)^2+(pi/b)^2),...
        fs(3) == c/(2*pi)*sqrt((pi/a)^2+(2*pi/b)^2)];
    s = solve(eqns);
end