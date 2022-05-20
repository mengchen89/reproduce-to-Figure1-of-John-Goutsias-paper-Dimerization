% parameters setting
A = 6.0221415*(10^23);
c1 = 10^-3;
V = 2;
%S = 10;
S = 1;
K = 6000;
s = S/(A*V);
k1 = A*V*c1;




t = linspace(1,7200,7200);

tend = 7200;
SS = [-1,-1,1];
x0 = [S,S,0];
Z_all = zeros(6000,7200);

%implement Gillespie algo.
for k = 1:K
    Z = [0];
    t_g = [0];
    while t_g(end)<tend
        current_z = Z(end);
        current_t = t_g(end);
        current_pp = a1_z(current_z,S);
        tau = exprnd(1/current_pp);
        Z_one = Z_all(k,:);
        Z_one(t >= current_t & t < current_t+tau) = current_z;
        Z_all(k,:) = Z_one;
        Z=[Z,Z(end)+1];
        t_g=[t_g,t_g(end)+tau];
    end
    Z_final = Z_all(k,:);
    Z_final(t >= current_t) = current_z;
    Z_all(k,:) = Z_final;
end

Z_ave = sum(Z_all)/K;

n = size(t);
I = n(2);
q3_array = [];
for i=1:I
    q3_array = [q3_array,q3_t(t(i),S)];
end


%derive population number
n_1 = size(Z_ave);
X1 = zeros(1,n_1(2));
X2 = zeros(1,n_1(2));

for i = 1:n_1(2)
    X1(i) = S-Z_ave(i);
    X2(i) = S-Z_ave(i);
end
X3 = Z_ave;
%derive mean concentration
u3 = X3/(A*V);

% plot of concentration
figure(1);
plot(t,u3);
hold on;
plot(t,q3_array);
xlabel('time'); 
ylabel('normalized dimer concentration');
legend({'SCKE','CKE'},'Location','southeast');

% derive mesoscopic forcing term and flux
v3 = std(Z_all,0,1).^2;
e3 = v3*k1/(A^2 * V^2);

n_2 = size(e3);
v1 = zeros(1,n_2(2));
for i = 1:n_2(2)
    v1(i) = k1*(s-u3(i))^2+e3(i);
end
% plot of mesoscopic forcing term
figure(2);
plot(t,e3);
xlabel('time'); 
ylabel('mesoscopic forcing term');

rhol_array = [];
for i=1:I
    rhol_array = [rhol_array,rhol_t(t(i),S)];
end
% plot of flux
figure(3)
plot(t,rhol_array);
hold on;
plot(t,v1);
xlabel('time'); 
ylabel('normalized flux');
legend({'SCKE','CKE'},'Location','northeast');

