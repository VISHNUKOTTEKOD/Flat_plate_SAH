function[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,X1,X2,X3,X4,X5]= coefff(T_g,T_air,T_p,T_b,T_ae,T_ag,dtau,dz,n,mdot,k) 
 
[delta_g,delta_p,delta_b,delta_air,delta_ag,c_g,c_p,c_b,c_ag,rho_g,rho_p,rho_b,rho_ag,alpha,tau_alpha,K_b,K_p,c_air]= get_constantsss;
[hcamg,hrpg,hcairp,hcagg,hcagp,hrpb,hcairb,Ub,hrbi]=get_h1(n,k,T_ae,T_g,T_p,T_ag,T_b,T_air,delta_b,K_b,dz,delta_ag);
[rho_a] = rhooo(T_air); 

A1=zeros(n,1);A2=zeros(n,1);A3=zeros(n,1);A4=zeros(n,1);A5=zeros(n,1);A6=zeros(n,1);A7=zeros(n,1);A8=zeros(n,1);A9=zeros(n,1);A10=zeros(n,1); 
A12=zeros(n,1);A13=zeros(n,1);A14=zeros(n,1);A15=zeros(n,1);A16=zeros(n,1);A17=zeros(n,1);X1=zeros(n,1);X2=zeros(n,1);X3=zeros(n,1);X4=zeros(n,1);X5=zeros(n,1);A11=zeros(n,1); 
 
for j=1:n

A1(j) = hcamg(j) / (rho_g * delta_g * c_g);

A2(j) = hcagg(j) / (rho_g * delta_g * c_g);
A3(j) = hrpg(j) / (rho_g * delta_g * c_g);
A16(j) = alpha / (rho_g * delta_g * c_g);
A4(j) = hcagg(j) / (rho_ag * delta_ag * c_ag);
A5(j) = hcagp(j) / (rho_ag * delta_ag * c_ag);
A6(j) = hrpg(j) / (rho_p * delta_p * c_p);
A7(j) = hcagp(j) / (rho_p * delta_p * c_p);
A8(j) = hcairp(j) / (rho_p * delta_p * c_p);
A9(j) = hrpb(j)/(rho_p * delta_p * c_p);
A17(j) = tau_alpha / (rho_p * delta_p * c_p);
A10(j) = hcairp(j) / (rhooo(T_air(j)) * delta_air * c_air);
A11(j) = hcairb(j) / (rhooo(T_air(j)) * delta_air * c_air);
A12(j) = mdot  / ( rhooo(T_air(j)) * delta_air  * 0.8*dz);
A13(j) = hrpb(j) / (rho_b * delta_b * c_b);
A14(j) = hcairb(j) / (rho_b * delta_b * c_b);
A15(j) = Ub(j) / (rho_b * delta_b * c_b);



X1(j) = (1/dtau) + A1(j) + A2(j) + A3(j)  ;
X2(j) = (1/dtau) + A4(j) + A5(j);
X3(j) = (1/dtau) + A6(j) + A7(j) + A8(j) + A9(j)  ;
X4(j) = (1/dtau) + A10(j) + A11(j) + A12(j)*dz;
X5(j) = (1/dtau) + A13(j) + A14(j) + A15(j);
end
