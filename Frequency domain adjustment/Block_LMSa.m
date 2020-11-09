clear;
n=153600; %%mikos tou u
M=1024;%%mikos tou filtrou
m=3.2e-04; %%o sintelestis m
L=M;%%epilogi tou block size isou me to megethos tou filtrou
k_max=n/M; %%arithmos twn block
%%upologismos tou u %%
v = sqrt(0.57)*randn(n,1); 
v = v - mean(v);  
u = zeros(n,1);
u(1) = v(1);
for i=2:n
  u(i) = -0.34 * u(i-1) + v(i);
end
%%%%%%%%%%%%%%%%%%%%%%%
 d=plant(u'); %%epithimito sima
 d=d';
 %%arxikopoiiseis%%
 y=zeros(n,1); %%eksodos tou filtrou
 e=zeros(n,1); %%sfalma
 w=zeros(L,1); %%filtro
 J=zeros(k_max,1); %%kampiles ekmathisis
 tic;
 for k=2:k_max %%epanalipsi gia kathe block
     %if k==2 %%tin prwti fora upologizetai to m
     %p=var(u((k-1)*L+1:k*L))*autocorr(u((k-1)*L+1:k*L),L-1);
     %m=0.9*(2/(L*eigs(toeplitz(p),1,'la')));
     %end
     phi=zeros(L,1);
     for i=1:L %%deuteri for loop
         for j=1:M
             y((k-1)*L+i)=y((k-1)*L+i)+w(j)*u((k-1)*L+i-j+1);%%upologismos tis eksodou tou filtrou
         end         
         e((k-1)*L+i)=d((k-1)*L+i)-y((k-1)*L+i); %%upologismos tou sfalmatos
         phi=phi+m*e((k-1)*L+i)*u((k-1)*L+i:-1:(k-2)*L+i+1); %%upologismos tou phi
         J(k)=J(k)+e((k-1)*L+i)^2;
     end
     w=w+phi; %%filtro
 end
 time=toc;
%% learning curves %%
figure
semilogy(J);
xlabel('Block Iterations')
ylabel('Error')
title('Learning curves, nested loops');