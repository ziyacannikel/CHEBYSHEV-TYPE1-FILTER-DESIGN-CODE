 clear all
 ws=0.5; % kesme bandý frekansý
 wp=0.7; % geçirme bandý frekansý
 rp=1; %geçirme bandý ripple
 rs=32; %kesme bandý ripple
 [b,a] = cheby1(n,Rp,Wp)
 % dijital verilen frekanslarý analog eksene çevirmek (prewarp) için
 wp=tan(pi*wp/2);
 ws=tan(pi*ws/2);

 
 w=wp/ws; % aþaðýdaki formülde karýþýklýk olmamasý için
 
 % chebyshev filtre derecesi bulmak için formül
 order=ceil(acosh(sqrt((10^(.1*rs)-1)/(10^(.1*rp)-1)))/acosh(w)); 
 
 
 wp=0.7; %yukarda wp bozuldu düzeltmek için
 
 %analog frekans elde etmek için formül
 fse = 2;
 u = 2*fse*tan(pi*wp/fse);
 

 %N. derece chebyshev analog prototip 
 [z,p,k] = cheb1ap(order,rp);
 
 %state-space forma geçiþ
 [a,b,c,d] = zp2ss(z,p,k);
 
 %alçak geçirenden yüksek geçiren filtreye geçiþ
 [a,b,c,d] = lp2hp(a,b,c,d,u);
 
 %biineer dönüþüm
 [a,b,c,d] = bilinear(a,b,c,d,fse);
 
 
 %----numaratör vektörünü kesinleþtirmek için bi yöntem---
 poles=poly(a);
 g=10^(-1/20); %% 1 = passband ripple
 Wn=0.7;
 Wn = 2*atan2(Wn,4);
 r=ones(4,1);%% matris oluþtur birlerden
 b=poly(r);
 kern = exp(-1i*pi*(0:length(b)-1));
 b = real(g*b*(kern*poles(:))/(kern*b(:)));
 zeros=b;
 %--------------------------------------------------------
 
 fs = 16e3; % örnekleme frekansý 
 dt = 1/fs; % örnekler arasý süre
 StopTime = 0.010; % grafik gösterim süresi
 t = (0:dt:StopTime)'; % süre 
 F = 5600; % Sine wave frequency (hertz) 
 data = sin(2*pi*F*t);% 6khz
 data2=sin(pi*F*t);% 3khz
 data3=data+data2;% sinüslerin toplamý
 
 subplot(2,2,1)
 plot(t,data,'-b');
 title('6khz sinüs'),xlabel('zaman'),ylabel('genlik')
 
 subplot(2,2,2)
 plot(t,data2,'-r');
 title('3khz sinüs'),xlabel('zaman'),ylabel('genlik')

 subplot(2,2,3)
 plot(t,data3)%filtresiz
 title('3khz+6khz sinüs'),xlabel('zaman'),ylabel('genlik')
 
 y=filter(zeros,poles,data3);
 
 subplot(2,2,4)
 plot(t,y,'-k');
 title('filtrelenmiþ sinyal'),xlabel('zaman'),ylabel('genlik')    
 



 
       
        
        
        
        
        
 


