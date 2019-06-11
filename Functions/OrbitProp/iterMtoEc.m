function Ec = iterMtoEc(M,e,tol)
er = 100;
Eco = M;
iter = 0;
while er >= tol && iter <= 100
    del = -(Eco-e*sin(Eco)-M)/(1-e*cos(Eco));
    Ec = Eco + del;
    er = abs(Ec-Eco);
    Eco = Ec;
    iter = iter + 1;
end