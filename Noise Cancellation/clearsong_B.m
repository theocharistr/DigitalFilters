clear all
close all

load('sound.mat') %%loading d
load('noise.mat') %%loading u

n=size(d,1);

%%pinakas autosusxetisi
n_aut = 3;
a = xcorr(u,u,n_aut-1,'unbiased');
a = a(n_aut:(2*n_aut-1));
R = toeplitz(a);

%%pinakas eterosusxetisis
P(1)=0.72;
P(2:3)=0;
P=P';
%%euresi idiotimwn kai m
idiot=eig(R);
max1=max(idiot);
%% fir filter

wo(1:3) =-1; 
w=wo; 
w=w'; %arxiki timi twn suntelestwn
mu = 2/max1;
mu=0.95*mu;

y = zeros(n, 1);

s = u;
for i=1:n
  w = w + mu*(P-R*w); % Adaptation steps
   if(i<n_aut)
    y(i)=s(i:-1:1)'*w(1:i);  %eksodos tou filtrou gia times tou i<3
    else
     y(i)=s(i:-1:i-n_aut+1)'*w; %filter
   end
end
song=d-y; %epithimiti eksodos 
