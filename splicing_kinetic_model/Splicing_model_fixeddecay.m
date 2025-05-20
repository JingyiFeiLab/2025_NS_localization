function y = Splicing_model_fixeddecay(t,para,var,deltap_fixed)
% This script implements a splicing kinetics model in which the pre-mRNA 
% splicing rate is fixed based on values obtained from previous experimental 
% measurements. 
ka=para(1);
deltap=deltap_fixed;
ks=para(2);
deltas=para(3);

beta1=ks+deltap;
P=(ka./beta1).*(1-exp((-beta1).*t));
S=(ka.*ks)./(beta1.*deltas)-(ka.*ks.*exp(-beta1.*t))./(beta1.*(deltas-beta1))+(ks.*ka.*exp(-deltas.*t))./(deltas.*(deltas-beta1));
y=[P;S];
y=y./var;
end