clear all;
close all;
%arxikopoiiseis metavlitwn
n=15000;
M=100;
nCoeff = 100;
s=zeros(n,1);
u=zeros(n,1);
x=zeros(n,1);
d=zeros(n,1);
y=zeros(n,1);
e=zeros(n,1);
J=zeros(n,1);
%upologismos simatwn
v1 = sqrt(0.42)*randn(n,1); 
v1 = v1 - mean(v1);
v2 = sqrt(0.72)*randn(n,1); 
v2 = v2 - mean(v2);
u(1:3)=v1(1:3);
for i=4:n
    u(i)=-0.87*u(i-1)-0.22*u(i-2)-0.032*u(i-3)+v1(i);
end
s(1:3)=u(1:3);
for i=4:n
    s(i)=-0.13*u(i)-0.67*u(i-1)-0.18*u(i-2)-0.39*u(i-3);
end
x(1:3)=v2(1:3);
for i=4:n
    x(i)=-0.57*x(i-1)-0.16*x(i-2)-0.08*x(i-3)+v2(i);
end
d(:)=s(:)+x(:);
%%%%%%%%%%%%%%%%%%%%%%%%
%upologismos pinaka autosisxetisis
[r_aut,lags] = xcorr(u,nCoeff-1,'unbiased');
r_aut=r_aut(lags>=0);
R = toeplitz(r_aut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upologismos pinaka eterosisxetisis
[p_et, lags] = xcorr(d, u, nCoeff-1, 'unbiased');
p_et = p_et(lags >= 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%veltistoi sintelestes wiener kai efarmogi
w1=zeros(M,1);
w1 =R\p_et;

y = zeros(n, 1);
e = zeros(n, 1);
J = zeros(n, 1);

for i=M+1:n
   y(i)=w1'*u(i:-1:i-M+1);
   e(i)=d(i)-y(i);
   J(i)=(e(i)-x(i))^2; %euresi sfalmatos
end
J_Wiener=J;
%RLS kai efarmogi
lambda=1;
P=250*eye(M,M);
wl = zeros(M, 1);
y = zeros(n, 1);
e = zeros(n, 1);
J = zeros(n, 1);
for i=M+1:n
    y(i)=wl(:)'*u(i:-1:i-M+1);
    k = ( (lambda^-1)*P*u(i:-1:i-M+1) / (1 + (lambda^-1)*u(i:-1:i-M+1)'*P*u(i:-1:(i-M+1))) );
    e(i) = d(i) - y(i);
    wl(:)=wl(:)+k*e(i);
    P = (lambda^-1)*P - (lambda^-1)*k*u(i:-1:(i-M+1))'*P;
    J(i) = (e(i)-x(i))^2;
end
J_RLS= mean(J,2);%upologismos sfalmatos

%LMS kai efarmogi
mu_LMS=0.00075;
wl = zeros(M, 1);
y = zeros(n, 1);
e = zeros(n, 1);
J = zeros(n, 1);

    for i = (M+1):n
        y(i) = wl(:)'*u(i:-1:(i-M+1));
        e(i) = d(i) - y(i);
        wl(:) = wl(:) + mu_LMS*e(i)*u(i:-1:i-M+1);
        J(i) = (e(i)-x(i))^2;
    end
J_LMS =J;%upologismos sfalmatos
%Normalized LMS kai efarmogi
wl = zeros(M, 1);
y = zeros(n, 1);
e = zeros(n, 1);
J = zeros(n, 1);

    for i = (M+1):n
        
    k = u(i:-1:(i-M+1))/ (e(i) + u(i:-1:i-M+1)'*u(i:-1:i-M+1) );
    e(i) = d(i) -wl(:)'*u(i:-1:(i-M+1)) ;
    wl(:)=wl(:)+k*e(i);
    
    J(i) = (e(i)-x(i))^2;
    end
 
J_Normalized=J;%upologismos sfalmatos
semilogy([J_LMS J_RLS J_Wiener J_Normalized])
legend({'LMS', 'RLS','Wiener','Norm'})
%sigkrisi LMS me Wiener
figure
semilogy([J_LMS J_Wiener])
legend({'LMS','Wiener'})

%sigkrisi RLS me Wiener
figure
semilogy([J_RLS J_Wiener])
legend({'RLS','Wiener'})

%sigkrisi Normalized LMS me Wiener
figure
semilogy([J_Normalized J_Wiener])
legend({'Normalized LMS','Wiener'})


