clear all
n=153600; %%mikos tou u
M=1024; %%mikos tou filtrou
mu=3.2e-04; %%o sintelestis m
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
alpha=0.5; %%step size
gamma=0.7; %%forgetting factor
P=1*ones(2*M,1); %%energy
w=zeros(2*M,1); %%filtro
d=d(:); %%epithimito sima
u=u(:); %%sima eisodou
e=d; %%sfalma
J=zeros(k_max,1); %%kampiles ekmathisis
tic;
for k=1:k_max-1 %%epanalipsi gia kathe block
 Uvector=fft([u((k-1)*M+1:(k+1)*M)],2*M); %%enwsi 2 block gia tin eisodo u kai euresi tou fft
 yvector=ifft(Uvector.*w); %%upologismos tis eksodou tou filtrou kai antistrofos fft
 yvector=yvector(M+1:2*M,1); %%kratietai mono to teleutaio block
 dvector=d(k*M+1:(k+1)*M); %%upologismos tou epithimitou simatos
 e(k*M+1:(k+1)*M,1)=dvector-yvector; %%upologismos sfalmatos
 Evector=fft([zeros(M,1);e(k*M+1:(k+1)*M)],2*M); %%euresi tou fft, gemisma twn prwtwn M stoixeiwn me 0
 P=gamma*P+(1-gamma)*abs(Uvector).^2; %%upologismos tou P
 Dvector=1./P;
 phi=ifft(Dvector.*conj(Uvector).*Evector,2*M); %%upologismos tou phi
 phi=phi(1:M);
 J(k)=J(k)+sum(real(dvector-yvector).^2); %%kampiles ekmathisis 
 w=w+alpha*fft([phi;zeros(M,1)],2*M); %%filtro
end
 e=real(e(:)); %%to sfalma prepei na exei mono pragmatikes times
 w=ifft(w); %%antistrofos fft gia epistrofi sto xrono
 w=real(w(1:length(w)/2));
 time=toc;

%% learning curves %%
figure
semilogy(J);
xlabel('Block Iterations')
ylabel('Error')
title('Learning curves,using FFT');