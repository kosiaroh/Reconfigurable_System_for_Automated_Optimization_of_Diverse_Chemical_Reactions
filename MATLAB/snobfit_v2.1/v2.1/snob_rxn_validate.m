function x = snob_rxn_validate(x,u,v)

n = length(x);
x(n/2+1:n) = sort(x(n/2+1:n),2,'descend');

if x(end)> v(end)
    x(end) = v(end);
elseif x(end) < u(end)
   x(end) = u(end);
end

for i = n-1:-1:n/2+1
    if x(i) < x(i+1) + 20
        x(i) = x(i+1) + 20;
    end
    if x(i) > v(i)
        x(i) = v(i);
    elseif x(i) < u(i)
        x(i) = u(i);
    end
end
end