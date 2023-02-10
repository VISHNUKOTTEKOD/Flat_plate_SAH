function[delta_g,delta_p,delta_b,delta_air,delta_ag,c_g,c_p,c_b,c_ag,rho_g,rho_p,rho_b,rho_ag,alpha,tau_alpha,K_b,K_p,c_air]= get_constantsss
delta_g=4/1000; % glass cover thickness(m)

delta_p=0.001; %plate thickness(m)
delta_b=0.01; %base thickness(m)
delta_air=.02; %air flow thickness(m)
delta_ag=0.04; %air gap thickness
c_g=800; %glass cover specific heat (J/kg.K)
c_ag=1000;
c_p=900; %plate specific heat (J/kg.K)
c_b=2300; %base insulation specific heat (J/kg.K)
rho_g=2500; %glass cover density(Kg/m^3)

rho_p=7850; %absorber density(Kg/m^3)
rho_b =22 ; %base plate density(Kg/m^3)
rho_ag =1.225 ; %assuume air gap density as constant 
alpha=.1; %absorption coefficient
tau_alpha=0.88; %trasmitivity*absorbtivity
K_b=0.14; %insulation thermal conductivity(W/m.K)
K_p=52; %absorber thermal conductivity(W/m.K)

c_air=1.0056e+003; %air specific heat (J/kg.K)

end