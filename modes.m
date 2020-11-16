function freqs = modes(c,a,b)
    freqs = [];
    Is = [];
    Js = [];
    for i = 1:100
        for j = 1:100
            freqs = [freqs, c * sqrt((i*pi/a)^2 + (j*pi/b)^2)];
            Is = [Is, i];
            Js = [Js, j];
        end
    end
    freqs = [freqs; Is; Js];
    freqs = freqs';
    freqs = sortrows(freqs);
    freqs(:,1) = freqs(:,1)/1000;
end