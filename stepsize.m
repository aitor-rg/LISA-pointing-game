function s = stepsize(t, R, a)
    if a<0.5
        a = 0.7;
    end
    s = 1/(t+R)^a;
end