clear all;
load speakerA.mat;
load speakerB.mat;
M=6600;
n=size(u,1);
w1=zeros(M,1);
wl=zeros(M,1);
y=zeros(n,1);
e=zeros(n,1);
%upologismos pinaka autosisxetisi
[r_aut,lags] = xcorr(u,M-1,'unbiased');
r_aut=r_aut(lags>=0);
R = toeplitz(r_aut);
%upologismos pinaka eterosisxetisis
[p_et, lags] = xcorr(d, u, M-1, 'unbiased');
p_et = p_et(lags >= 0);
%veltistoi sintelestes wiener kai efarmogi
w1 =R\p_et;
%me veltistous wiener sintelestes 
for i=M+1:n
   y(i)=w1'*u(i:-1:i-M+1);
   e(i)=d(i)-y(i);
end
%USING LMS ALGORITHM 
%{
mu=0.001;

for i = M+1:n
        y(i) = wl(:)'*u(i:-1:(i-M+1));
        e(i) = d(i) - y(i);
        wl(:) = wl(:) + mu*e(i)*u(i:-1:i-M+1);        
        if (isnan(y(i)) || isnan(e(i)))
            break;
        end
end
%}

