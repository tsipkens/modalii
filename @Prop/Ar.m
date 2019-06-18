function [] = Ar(prop)

prop.mg = 6.63368e-26;
prop.gamma1 = 5/3;
prop.gamma2 = @(T)(prop.gamma1+1)/(prop.gamma1-1);

end

