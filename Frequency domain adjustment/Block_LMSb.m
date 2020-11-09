clear all
n=153600; %%mikos tou u
M=1024; %%mikos tou filtrou
k_max=n/M; %%arithmos twn block
mu=3.2e-04; %%o sintelestis m
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
w=zeros(M,1); %%filtro
J=zeros(k_max,1); %%kampiles ekmathisis
tic;
for k=1:k_max-1 %%epanalipsi gia kathe block
    %if k==1 %%tin prwti fora upologizetai to m
    %p=var(u((k)*M+1:(k+1)*M))*autocorr(u((k)*M+1:(k+1)*M),M-1);
    %mu=0.7*(2/(M*eigs(toeplitz(p),1,'la')));
    %end
umatrix=toeplitz(u(k*M:1:(k+1)*M-1),u(k*M:-1:(k-1)*M+1)); %% euresi tou pinaka pou periexei tin eisodo
dvector=d(k*M:1:(k+1)*M-1); %%euresi tou dianismatos pou periexei to epithimito sima
yvec=umatrix*w; %%upologismos simatos eksodou
evector=dvector-yvec; %%upologismos dianismatos sfalmatos
e(k*M:1:(k+1)*M-1)=evector; %%apothikeusi ston katallilo pinaka e
phi=umatrix.'*evector; %%upologismos tou phi
w=w+mu*phi; %%filtro
J(k)=J(k)+sum(evector.^2);  
end
time=toc;
%% learning curves %%
figure
semilogy(J);
xlabel('Block Iterations')
ylabel('Error')
title('Learning curves,one loop');