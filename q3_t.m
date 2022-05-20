function q3t = q3_t(t,S)
A = 6.0221415*(10^23);
c1 = 10^-3;
V = 2;
K = 6000;
s = S/(A*V);
k1 = A*V*c1;
q3t = k1*(s^2)*t/(1+k1*s*t); 
end