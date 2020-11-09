function y = T(n)
if n==2
    y=4; %%logw tou oti T(2)=4
else
    y=5*n+2*T(n/2); %%anadromiki sxesi
end
end

