function k_a=ka(T_air)
k_a=((0.0015215)+(0.097459.*T_air)-((3.3322*10^(-5))*(T_air.^2)))/(1000);
end