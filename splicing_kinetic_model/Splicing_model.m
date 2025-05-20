function y = Splicing_model(t,para,var)
% This script implements the full splicing kinetics model where both the 
% splicing rate and degradation rate are treated as free parameters to be 
% fitted from the data. 
% This model can be directly compared to the fixed-rate version 
% (Splicing_model_fixeddecay) 
ka=para(1);
deltap=para(2);
ks=para(3);
deltas=para(4);

beta1=ks+deltap;
P=(ka./beta1).*(1-exp((-beta1).*t));
S=(ka.*ks)./(beta1.*deltas)-(ka.*ks.*exp(-beta1.*t))./(beta1.*(deltas-beta1))+(ks.*ka.*exp(-deltas.*t))./(deltas.*(deltas-beta1));
y=[P;S];
y=y./var;
end