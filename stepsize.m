function s = stepsize(t, R, a,g)
%     if a<0.5
%         a = 0.7;
%     end
    s = g/(t+R)^a;
end