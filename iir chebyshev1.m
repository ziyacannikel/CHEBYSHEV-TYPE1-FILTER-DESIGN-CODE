 clear all
 ws=0.5; % kesme band� frekans�
 wp=0.7; % ge�irme band� frekans�
 rp=1; %ge�irme band� ripple
 rs=32; %kesme band� ripple
 [b,a] = cheby1(n,Rp,Wp)
 % dijital verilen frekanslar� analog eksene �evirmek (prewarp) i�in
 wp=tan(pi*wp/2);
 ws=tan(pi*ws/2);

 
 w=wp/ws; % a�a��daki form�lde kar���kl�k olmamas� i�in
 
 % chebyshev filtre derecesi bulmak i�in form�l
 order=ceil(acosh(sqrt((10^(.1*rs)-1)/(10^(.1*rp)-1)))/acosh(w)); 
 
 
 wp=0.7; %yukarda wp bozuldu d�zeltmek i�in
 
 %analog frekans elde etmek i�in form�l
 fse = 2;
 u = 2*fse*tan(pi*wp/fse);
 

 %N. derece chebyshev analog prototip 
 [z,p,k] = cheb1ap(order,rp);
 
 %state-space forma ge�i�
 [a,b,c,d] = zp2ss(z,p,k);
 
 %al�ak ge�irenden y�ksek ge�iren filtreye ge�i�
 [a,b,c,d] = lp2hp(a,b,c,d,u);
 
 %biineer d�n���m
 [a,b,c,d] = bilinear(a,b,c,d,fse);
 
 
 %----numarat�r vekt�r�n� kesinle�tirmek i�in bi y�ntem---
 poles=poly(a);
 g=10^(-1/20); %% 1 = passband ripple
 Wn=0.7;
 Wn = 2*atan2(Wn,4);
 r=ones(4,1);%% matris olu�tur birlerden
 b=poly(r);
 kern = exp(-1i*pi*(0:length(b)-1));
 b = real(g*b*(kern*poles(:))/(kern*b(:)));
 zeros=b;
 %--------------------------------------------------------
 
 fs = 16e3; % �rnekleme frekans� 
 dt = 1/fs; % �rnekler aras� s�re
 StopTime = 0.010; % grafik g�sterim s�resi
 t = (0:dt:StopTime)'; % s�re 
 F = 5600; % Sine wave frequency (hertz) 
 data = sin(2*pi*F*t);% 6khz
 data2=sin(pi*F*t);% 3khz
 data3=data+data2;% sin�slerin toplam�
 
 subplot(2,2,1)
 plot(t,data,'-b');
 title('6khz sin�s'),xlabel('zaman'),ylabel('genlik')
 
 subplot(2,2,2)
 plot(t,data2,'-r');
 title('3khz sin�s'),xlabel('zaman'),ylabel('genlik')

 subplot(2,2,3)
 plot(t,data3)%filtresiz
 title('3khz+6khz sin�s'),xlabel('zaman'),ylabel('genlik')
 
 y=filter(zeros,poles,data3);
 
 subplot(2,2,4)
 plot(t,y,'-k');
 title('filtrelenmi� sinyal'),xlabel('zaman'),ylabel('genlik')    
 



 
       
        
        
        
        
        
 


