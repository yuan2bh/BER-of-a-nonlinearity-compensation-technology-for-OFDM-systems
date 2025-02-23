function [x_GI,OBO] = nonlPA(x_GI_temp,dispTemp,amp,bq_CO,P)

Lbq = length(bq_CO);   

% determine the OBO of PA
% amp = 20;

dn = x_GI_temp * amp;
Ldn=length(dn);

[Dson,GMPmode_D] = P_onlyOdd_Dson(dn,P);

xn=Dson*bq_CO;     %PA output
xnT=xn.';

OBO=10*log10( max(abs(xn).^2) / mean(abs(xn).^2));


if(dispTemp)
figure;
subplot(221);    plot(abs(x_GI_temp),abs(xnT),'*');  title('AM-AM');     hold on;
ang = angle(x_GI_temp) - angle(xnT);
subplot(223);    plot( abs(x_GI_temp), unwrap(ang),'*');  title('AM-PM');     hold on;
end  
x_GI = xnT;



