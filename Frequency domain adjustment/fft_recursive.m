function Y = fft_recursive(x)
n = length(x); %%mikos tis eisodou 
global muls; %%o arithmos twn pol/smwn
global sums; %%o arithmos twn prosthesewn
if (n == 1)
  Y = conj(x);
else 
  Y1 = fft_recursive(x(1:2:n)); %%peritta stoixeia
  m = n/2;
  W = exp(2*pi*1i*[0:1:m-1]/n);
  Y2 = fft_recursive(x(2:2:n)).*W; %%zuga stoixeia
  %% upologismos twn prosthesewn kai pol/smwn%%
  if n==2
      sums=sums+2; %%T(2)=4
  else
  muls=muls+n/2;
  sums=sums+2*length(Y2);
  end
  %% TELOS %%
  Y = [Y1+Y2 Y1-Y2];
end