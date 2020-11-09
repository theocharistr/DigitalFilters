close all;
clear 

%arxikopoiiseis metablitwn%
n=1000;
nCoeff = 100;
A=4.2;
fo=1/4;
phi=pi/2;
Delta=10;
J=zeros(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upologismos simatwn
v1 = sqrt(0.54)*randn(n,1); 
v1 = v1 - mean(v1);
for i=1:n
x(i)=A*(sin(2*pi*fo*i+phi)+cos(4*pi*fo*i+phi)+cos(7*pi*i+phi/3));
end
s=v1+x';
u=zeros(n,1);
for i=Delta+1:n
    u(i)=s(i-Delta);
end
d=s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%upologismos pinaka autosisxetisis
[r_aut,lags] = xcorr(u,nCoeff-1,'unbiased');
r_aut=r_aut(lags>=0);
R = toeplitz(r_aut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upologismos pinaka eterosisxetisis
[p_et, lags] = xcorr(s, u, nCoeff-1, 'unbiased');
p_et = p_et(lags >= 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%veltistoi sintelestes wiener kai efarmogi
e=zeros(n,1);
w1 =R\p_et;
y(1:nCoeff)=u(1:nCoeff);
for i=nCoeff+1:n
   y(i)=w1'*u(i:-1:i-nCoeff+1);
   e(i)=s(i)-y(i);
   J(i)=(e(i)-v1(i))^2; %upologismos tou sfalmatos
end
J_Wiener=J;
%klisi tou algorithmou Levinson Durvin
[a, G, L, Dp]=LevinsonDurbin_iterative(nCoeff-1, r_aut);
[a1,e1,k1]=levinson(r_aut,nCoeff-1);
a1=a1';
k1=k1';
for i=1:nCoeff
    sfalma(i)=(a1(i)-a(i))^2;
end
sfalma_1=mean(sfalma); %sfalma metaksi twn dio klisewn 
%ektelesi tou lattice
b=zeros(n,nCoeff);
f=zeros(n,nCoeff);
b(:,1)=u(:,1);
f(:,1)=u(:,1);
for i=2:nCoeff
    for j=2:n
        f(j,i)=f(j,i-1)+G(i-1)*b(j-1,i-1);
        b(j,i)=G(i-1)*f(j,i-1)+b(j-1,i-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%ektelesi lattice me etoimi sinartisi tou matlab
[f_opt,b_opt] = latcfilt(G,u);
for i=1:n
    sfalma_f(i)=(f_opt(i)-f(i,100))^2;
end
sfalma_2=mean(sfalma_f); %sfalma metaksi twn dio klisewn 

D=zeros(nCoeff);
%%%%%%%%%%%%%%%%%%%%%%%%
%upologismos tou D gia tin euresi tou gamma
for i=1:nCoeff
    for j=1:nCoeff
        if i==j
            D(i,i)=Dp(i);
        end 
    end
end
gamma=inv(D)*L*p_et;%euresi mesw tou tipou
gamma_opt=L'\w1;%euresi me xrisi twn veltistwn sintelestwn wiener
for i=1:nCoeff
    sfalma_gamma(i)=(gamma_opt(i)-gamma(i))^2;
end
sfalma_3=mean(sfalma_gamma); %sfalma metaksi twn sintelestwn gamma
%%%%%%%%%%%%%%%%%%%%%%%%
%joint process estimation
y=zeros(n,1);
e=zeros(n,1);
J=zeros(n,1);
for i=1:n
    for j=1:nCoeff
        y(i)=y(i)+b(i,j)*gamma(j);
    end
    e(i)=s(i)-y(i);
    J(i)=(e(i)-v1(i))^2;
end
J_Joint=J;
figure
semilogy(J_Wiener );
legend('Optimal Wiener');
figure
semilogy(J_Joint);
legend('Joint process estimation');

