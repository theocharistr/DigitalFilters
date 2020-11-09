clear all;
%arxikopoiiseis kai upologismos tou u se sxesi me to s
load music.mat
nCoeff=100;
n=size(s,1);
u=zeros(n,1);
D=100;
for i=1:n
    if i<D+1
        u(i)=s(i);
    else
        u(i)=s(i-D);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upologismos pinaka autosisxetisis
[r_aut,lags] = xcorr(u,u,nCoeff-1,'unbiased');
r_aut = r_aut(lags>=0);
R = toeplitz(r_aut);
%upologismos pinaka eterosisxetisis
[p_et,lags] = xcorr(u,s,nCoeff-1,'unbiased');
p_et=p_et(lags>=0);
%veltistoi sintelestes wiener kai efarmogi tous
w=R\p_et;
%{
y(1:nCoeff)=u(1:nCoeff);
for i=nCoeff+1:n
   y(i)=w'*u(i:-1:i-nCoeff+1);
   e(i)=s(i)-y(i);
  
end
%}
%klisi tou algorithmou Levinson Durvin
[a, G, L, Dp]=LevinsonDurbin_iterative(nCoeff-1, r_aut);
%lattice n=4000000 poli megalo dhmiourgei problhma
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
%%%%%%%%%%%%%%%%%%%%%%
%joint process estimator
gamma=L'\w;
y=zeros(n,1);
e=zeros(n,1);
J=zeros(n,1);
for i=1:n
    for j=1:nCoeff
        y(i)=y(i)+b(i,j)*gamma(j);
    end
    e(i)=s(i)-y(i);
    J(i)=e(i)^2;
end

