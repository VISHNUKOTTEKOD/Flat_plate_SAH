function Thermal_Air_Diffusivity= tad(Tm)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Thermal_Air_Diffusivity=((0.0146*(Tm-273.15))+(1.8343))/(10^4);
end