function D = objectD(L, w1, w2)
    f = 200;
    imgD = L/(1-w2/w1);
    D = 1/(1/f-1/imgD);
end