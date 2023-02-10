
ts=cputime;
n=2;
interval=720;
L=1;
w=1; delta_air=.02;
mdot=0.025;
%V=mdot/1.225;
%v=V/(w*delta_air);
v=1;
dz=L/(n-1);
dtau=1;
T_int=interval*60;
T_tot=round(T_int/dtau);

if dtau>dz/v
   fprintf('error in flow rate')
else 
    tfile = 'ptemp.out';
    fid = fopen(tfile,'wt');
    T_ae=zeros(T_tot+1,1);
    G_r=zeros(T_tot+1,1);
    T_g=ones(n,1)*294.66;
    %T_s=ones(n,1)*294.66;
    T_p=ones(n,1)*294.66;
    T_air=ones(n,1)*294.66;
    T_b=ones(n,1)*294.66;
    T_ag=ones(n,1)*294.66;
    T_gc=zeros(T_tot,1);T_agc=zeros(T_tot,1);
    T_pc=zeros(T_tot,1);T_airc=zeros(T_tot,1);
    T_bc=zeros(T_tot,1);T_out=zeros(T_tot,1);
    Q_dot=zeros(T_tot,1);
    Q_dot2=zeros(T_tot,1);
    eff=zeros(T_tot,1);
    peff=zeros(T_tot/2,1);
    
   

    for k=1:T_tot   %%ttot to ttot+1
        if k*dtau<=43200    %%%
           G_r(k)=950*sin((3.14*k)/43200);
          %T_ae(k)=(-0.00009*(k/60)^2)+(0.0743*k/60)+294.66;
          
            %%%%%%%%%%%%%%
            openfig('950c.fig')
            ax = gca; 
            h = findobj(gca,'Type','line');

            k = h.XData;
            y = h.YData;

            %
            P = polyfit(k,y,2); % 2 is the degree of polynomial fit.If you increase that, the accuracy of fit increases
             %but the simpler the best for analysis

            y1 = polyval(P,k); % to draw the fitted y values
            T_ae(k)=y1(k)';
            %%%%%%%%%%%%%%%%
          
           %T_ae(k)=(-5*(10^-5)*(k/60)^2+0.0468*(k/60)+294.66);
        else
        G_r(k)=0;
        T_ae(k)=(-0.00009*(k/60)^2)+(0.0743*k/60)+294.66;
        end
    end
    for k=1:T_tot
        n_converge=0;
       [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,X1,X2,X3,X4,X5]= coefff(T_g,T_air,T_p,T_b,T_ae,T_ag,dtau,dz,n,mdot,k);
       kk=0;
       while n_converge<5*n
           kk=kk+1;
           T_g_old=T_g;T_ag_old=T_ag;T_p_old=T_p;T_air_old=T_air;T_b_old=T_b;
           
           T_g(1)=((T_g_old(1)/dtau)+(A1(1)*T_ae(k))+(A2(1)*T_ag(1))+(A3(1)*T_p(1))+(A16(1)*G_r(k)))/X1(1);
           T_ag(1)=((T_ag_old(1)/dtau)+(A4(1)*T_g(1))+(A5(1)*T_p(1)))/X2(1);
           T_p(1)=((T_p_old(1)/dtau)+(A6(1)*T_g(1))+(A7(1)*T_ag(1))+(A8(1)*T_air(1))+(A9(1)*T_b(1)+(A17(1)*G_r(k))))/X3(1);
           T_air(1)=T_ae(k);
           T_b(1)=((T_b_old(1)/dtau)+(A13(1)*T_p(1))+(A14(1)*T_air(1))+(A15(1)*T_ae(k)))/X5(1);

           for j=2:n 
            T_g(j)=((T_g_old(j)/dtau)+(A1(j)*T_ae(k))+(A2(j)*T_ag(j))+(A3(j)*T_p(j))+(A16(j)*G_r(k)))/X1(j);
           T_ag(j)=((T_ag_old(j)/dtau)+(A4(j)*T_g(j))+(A5(j)*T_p(j)))/X2(j);
           T_p(j)=((T_p_old(j)/dtau)+(A6(j)*T_g(j))+(A7(j)*T_ag(j))+(A8(j)*T_air(j))+(A9(j)*T_b(j)+(A17(j)*G_r(k))))/X3(j);
           T_air(j)=((T_air_old(j)/dtau)+(A10(j)*T_p(j))+((A11(j)/dz)*T_air(j-1))+(A12(j)*T_b(j)))/X4(j);
           T_b(j)=((T_b_old(j)/dtau)+(A13(j)*T_p(j))+(A14(j)*T_air(j))+(A15(j)*T_ae(k)))/X5(j);  
           end   
            %T_g(n)=((T_g_old(n)/dtau)+(A1(n)*T_ae(k))+(A2(n)*T_p(n))+(A3(n)*T_air(n))+(A4(n)*G_r(k)))/X1(n);
           %T_air(n)=((T_air_old(n)/dtau)+(A5(n)*T_g(n))+(A6(n)*T_p(n))+((A7(n)/dz)*T_air(n-1)))/X2(n);
           %T_p(n)=((T_p_old(n)/dtau)+(A8(n)*T_air(n))+(A9(n)*T_g(n))+(A10(n)*T_s(n))+(A11(n)*G_r(k)))/X3(n);
           %T_s(n)=((T_s_old(n)/dtau)+(A12(n)*T_p(n))+(A13(n)*T_b(n)))/X4(n);
           %T_b(n)=((T_b_old(n)/dtau)+(A14(n)*T_s(n))+(A15(n)*T_ae(k)))/X5(n);  
           

    %        end
           %check convergence
           ccc=0;
           for j=1:n
               if ccc<=0
                   error=zeros(5,1);
                   error(1)=abs(T_g(j)-T_g_old(j))/T_g(j);
                   error(2)=abs(T_ag(j)-T_ag_old(j))/T_ag(j);
                   error(3)=abs(T_p(j)-T_p_old(j))/T_p(j);
                   error(4)=abs(T_air(j)-T_air_old(j))/T_air(j);
                   error(5)=abs(T_b(j)-T_b_old(j))/T_b(j); 
                   for i=1:5
                       if(error(i)<=10^-4)
                           n_converge=n_converge+1;
                       else
                           ccc=1;
                       end
                   end
               end
           end 
       end
       T_gc(k)=T_g(n/2);T_agc(k)=T_ag(n/2);%T_airc(k)=T_air(n/2);
       T_pc(k)=T_p(n/2);T_bc(k)=T_b(n/2);
       T_out(k)=T_air(n);


       Q_dot(k)= 00.025*1005.6*(T_out(k)-T_ae(k));
      eff(k)=(Q_dot(k)/(1*0.8*G_r(k)))*100;

        time = dtau*k/60;
       fprintf('time = %6.1f T_out = %8.2f T_gu = %8.2f T_p = %8.2f T_ag = %8.2f T_b = %8.2f  n_iter =%5.0f\n',time,T_out(k),T_gc(k),T_pc(k),T_agc(k),T_bc(k),kk);
       fprintf('A1 = %6.8f A2 = %8.8f A3 = %8.8f A16 = %8.8f A4 = %8.8f A5 = %8.8f A6 = %8.8f A7 = %8.8f A8 =%8.8f A9 = %8.8f A17 = %8.8f, A10 = %8.8f A11 = %8.8f A12 = %8.8f A13 = %8.8f A14 = %8.8f A15 = %8.8f X1 = %8.8f X2 = %8.8f X3 = %8.2f X4 = %8.2f X5 = %8.2f \n',A1(j),A2(j),A3(j),A16(j),A4(j),A5(j),A6(j),A7(j),A8(j),A9(j),A17(j),A10(j),A11(j),A12(j),A13(j),A14(j),A15(j),X1(j),X2(j),X3(j),X4(j),X5(j));
        count=fprintf(fid,'%6.1f %6.1f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',time,G_r(k),T_out(k),T_gc(k),T_airc(k),T_pc(k),T_agc(k),T_bc(k));

    end
    for k=1:T_tot/2
        if G_r(k)<=100
            peff(k)=0;
        else
            peff(k)=eff(k);
        end
    end    

    T=6:12/43199:18;
    T2=1:(T_tot/2);
    le=1:n;
    figure(1)
    plot(T,T_gc,'r',T,T_pc,'y',T,T_agc,'g',T,T_bc,'k',T,T_ae(1:T_tot),'--')
    xlabel('Time(h)')
    ylabel('Temperature(K)')
    legend('T gc','T pc','T agc','T bc','T ae')

    figure(2)
    plot(T,T_ae(1:T_tot),'--',T,T_out)
    xlabel('Time(s)')
    ylabel('Tout(K)')
    legend('T ae','T out')
    figure(3)
    plot(T2,peff)
    xlabel('Time(s)')
    ylabel('Efficiency(%)')
    runtime = cputime-ts;
    
    figure(4)
    plot(T,G_r(1:T_tot))

    status = fclose(fid);
    Values=[max(T_gc);max(T_airc);max(T_pc);max(T_agc);max(T_bc);max(T_out);max(T_ae);peff(k)];
end




