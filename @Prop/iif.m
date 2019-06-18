function out = iif(cond,a,b)
    out = b;
    out(cond) = a(cond);
end