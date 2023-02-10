function[hcamg,hrpg,hcairp,hcagg,hcagp,hrpb,hcairb,Ub,hrbi]=get_h1(n,k,T_ae,T_g,T_p,T_ag,T_b,T_air,delta_b,K_b,dz,delta_ag)
hcamg=zeros(n,1);hrpg=zeros(n,1);hcairp=zeros(n,1);
hcagp=zeros(n,1);hcagg=zeros(n,1);hrpb=zeros(n,1);
%hrpb=zeros(n,1);hcair2p=zeros(n,1);hcair2b=zeros(n,1);
hcairb=zeros(n,1);Ub=zeros(n,1);
hrbi=zeros(n+1,1);
Ra=zeros(n,1);Nu_a=zeros(n,1);
Re2=zeros(n,1);Nu2=zeros(n,1);
Re1=zeros(n,1);Nu1=zeros(n,1);

%Re3=zeros(n,1);Nu3=zeros(n,1);
%Re4=zeros(n,1);Nu4=zeros(n,1);
g=9.81;w=0.8;th1=0.05;L=1.8;mdot=0.03;%th2=0.05;
vext=1.3;
T_sky(k)=0.0522*(T_ae(k).^1.5);
hw=5.67+(3.86*vext);
segma=5.667*10^-8;
theta=pi*11/180;
Dh1=(2*w*th1)/(w+th1);
%Dh2=(2*w*th2)/(w+th2);
emi_b=0.93;
emi_g=0.89;
emi_p=0.9;
%emi_b=0.9;
for j=1:n
     Ra(j)=(abs(T_g(j)-T_p(j))*g*L^3*rhooo(T_ag(j)))/(dv(T_ag(j))*tad(T_ag(j))*T_ag(j));  %%thickness to length
        AA=1-(1708*(sin(1.8*theta))^1.6/(Ra(j)*cos(theta)));
        BB=(Ra(j)*cos(theta)/5830)^(1/3)-1;
        if AA<=0
            if BB<=0
                Nu_a(j)=1;
            else 
                Nu_a(j)=1+BB;
            end
        else 
            if BB<=0
            Nu_a(j)=1+(1.44*(1-(1708*(sin(1.8*theta))^1.6/(Ra(j)*cos(theta))))*AA);
            else
            Nu_a(j)=1+(1.44*(1-(1708*(sin(1.8*theta))^1.6/(Ra(j)*cos(theta))))*AA)+BB;
            end
        end
    hcagg(j)=Nu_a(j)*ka(T_ag(j))/delta_ag;
    hcagp(j)=hcagg(j);   
    if T_g(j)-T_ae(k)==0
        hcamg(j)=hw;
        else
        hcamg(j)=((segma*emi_g*(T_g(j).^4-T_sky(k).^4))./(T_g(j)-T_sky(k)))+hw;
    end
    hrpg(j)=(segma*(T_g(j).^2+T_p(j).^2)*(T_g(j)+T_p(j)))/((1/emi_g)+(1/emi_p)-1);
    %Hc of air gap(pending)--hcaggu hcaggl
    hrpb(j)=(segma*(T_p(j).^2+T_b(j).^2)*(T_p(j)+T_b(j)))/((1/emi_p)+(1/emi_b)-1);
    %hrpb(j)=(segma*(T_p(j).^2+T_b(j).^2)*(T_p(j)+T_b(j)))/((1/emi_b)+(1/emi_b)-1);
    Tm1=(T_p(j)+T_air(j))/2;  %%Tg to Tp
    Re1(j)=(mdot*Dh1)/(w*th1*dv(Tm1)); %%m to mdot
    if(Re1(j)>2100)
    Nu1(j)=0.0158*(Re1(j)^0.8);   
    end
    if(Re1(j)<=2100)
        Num1=0.0606*((Re1(j)*Pr(Tm1)*Dh1/dz)^1.2);
        Den1=0.0909*((Re1(j)*Pr(Tm1)*Dh1)^0.7/(dz^1.2))*(Pr(Tm1)^0.7);
        Nu1(j)=4.9+((Num1)/(1+Den1));
    end
    if(Re1(j)>2100)
        Nu1(j)=0.0158*Re1(j)^0.8;
    end
    hcairp(j)=((Nu1(j)*ka(Tm1))/Dh1);
    Tm2=(T_b(j)+T_air(j))/2;  %%Tp to Tb
    Re2(j)=(mdot*Dh1)/(w*th1*dv(Tm2));
    if(Re2(j)>2100)
        Nu2(j)=0.0158*(Re2(j)^0.8);
    end
    if(Re2(j)<=2100)
        Num2=0.0606*((Re2(j)*Pr(Tm2)*Dh1/dz)^1.2);
        Den2=0.0909*((Re2(j)*Pr(Tm2)*Dh1)^0.7/(dz^1.2))*(Pr(Tm2)^0.7);
        Nu2(j)=4.9+((Num2)/(1+Den2));
    end
    hcairb(j)=((Nu2(j)*ka(Tm2))/Dh1);
        %Tm3=(T_air2(j)+T_p(j))/2;
        %Re3(j)=(m*Dh2)/(w*th2*dv(Tm3));
        %if(Re3(j)>2100)
        %Nu3(j)=0.0158*(Re3(j)^0.8);   
        %end
        %if(Re3(j)<=2100)
        %Num3=0.0606*((Re3(j)*Pr(Tm3)*Dh2/dz)^1.2);
        %Den3=0.0909*((Re3(j)*Pr(Tm3)*Dh2/dz)^0.7)*(Pr(Tm3)^0.17);
        %Nu3(j)=4.9+((Num3)/(1+Den3));
        %end
        %hcair2p(j)=((Nu3(j)*ka(Tm3))/Dh2);
        %Tm4=(T_air2(j)+T_b(j))/2;
        %Re4(j)=(m*Dh2)/(w*th2*dv(Tm4));
        %if(Re4(j)>2100)
        %Nu4(j)=0.0158*(Re4(j)^0.8);   
        %end
        %if(Re4(j)<=2100)
        %Num4=0.0606*((Re4(j)*Pr(Tm4)*Dh2/dz)^1.2);
        %Den4=0.0909*((Re4(j)*Pr(Tm4)*Dh2/dz)^0.7)*(Pr(Tm4)^0.17);
        %Nu4(j)=4.9+((Num4)/(1+Den4));
        %end
        %hcair2b(j)=((Nu4(j)*ka(Tm4))/Dh2);
    if T_b(j)-T_ae(k)==0
            hrbi(j)=hw;
    else
            hrbi(j)=((segma*emi_b*(T_b(j)^4-T_sky(k)^4))/(T_b(j)-T_sky(k)))+hw;
    end
       Ub(j)=(1)/((delta_b/K_b)+(1/(hw+hrbi(j))));
end
end





