function rholt = rhol_t(t,S)
A = 6.0221415*(10^23);
c1 = 10^-3;
V = 2;
K = 6000;
s = S/(A*V);
k1 = A*V*c1;
rholt = k1*(s^2)/(1+k1*s*t)^2; 
end