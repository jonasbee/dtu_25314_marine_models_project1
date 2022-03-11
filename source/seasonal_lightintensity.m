function I = seasonal_lightintensity(P,param,t)

% integral
int = cumsum(P.*param.dz)*param.k;

I = param.I_0*exp(-param.K_bg*param.z-int)*(1-cos(2*pi*t/365));

end