clear all
global sums;
global muls;
sums=0;
muls=0;
x =[-0.1241 + 0.4889i;
   1.4897 + 1.0347i;
   1.4090 + 0.7269i;
   1.4172 - 0.3034i;
   0.6715 + 0.2939i;
  -1.2075 - 0.7873i;
   0.7172 + 0.8884i;
   1.6302 - 1.1471i];
y=fft_recursive(x);
compl=2*sums+6*muls; %%poliplokotita mesw tou anadromikou upologismou tou fft
compl_2=T(length(x)); %%poliplokotita mesw tis anadromikis sxesis pou dinetai
ys=fft(x); %%ypologismos fft 
e=norm(y-ys'); %% euresi sfalmatos
if compl==compl_2
    fprintf('result confirmed with recursive fft algorithm error =: %e\n', e)
end