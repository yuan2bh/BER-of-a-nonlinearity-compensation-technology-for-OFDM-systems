% run ok @ MATLAB R2013b. 
% Honglin Yuan (Ԭ����)
% Nantong University, P.R.China (�й���ͨ��ѧ)
% 2025.2.20
%
%
% �����Ӽ�����BERͼ����2015.11.28.
%

clear all;  clc;
warning('off','comm:obsolete:randint');

% 
Nonliear = 1

% �����źŹ��ʿ���
meaSwitch = 0

% Number of iterations for each EbN0
N_iter=1e5;       

% Number of OFDM symbols per frame
% һframe�ڣ��� Nframe=Nesti+Ncomp��OFDM����
% ��Ϊ�����֣�ǰ��Nesti�����ڷ������������Ncomp�����ڲ���& BERͳ��
% ����Nesti��OFDM���Ž��з�����ϵ���������ŵ�������Ӧ�ķ��룻
% ����Ncomp�����Ž��л��ڷ�������Ĳ�����BERͳ��
Nframe= 103         
Nesti = 100     % ���ڲ������� 

% Ch=0/1 for AWGN/multipath channel���ŵ����ͣ�AWGN��ྶ��
Ch= 1  
if Ch==0
    chType='AWGN'; 
    Target_neb=100; 
else
    chType='CH'; 
    Target_neb=500*10*2*2; 
end

% PA�Ļ����Ŵ�����Ӱ��sigPow
PAbasicAmp=15

if(meaSwitch==1)
    sigPow = 0;
else
    if(Ch==0)
        sigPow =2.762e-03           % only ������ �²������ + Nframe=3
    else
        sigPow = 2.3e-1             % PAbasicAmp=15  &&   Nframe= 103 
    end
end

EbN0=[0:5:25];    % EbN0

% NgType=1/2 for cyclic prefix��ѭ��ǰ׺��/zero padding   
NgType=1 
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end

% ********************* chanenl ******************************
PowerdB=[0 -8 -17]; % Channel tap power profile 'dB'
Delay=[0 1 2];          % Channel delay 'sample'

Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Ntap=length(PowerdB);       % Chanel tap number
Lch=Delay(end)+1;           %Channel length

% PA
bq_CO=[1; -0.2+0.01*i; 0.1+0.05*i; 0.12-0.2*i ];
P=7;
Lbq=length(bq_CO);

% ********************** sigPow=1.39e-2 *********************************
Nbps=4; 
Nfft=1024;           % FFT size

% ************************************************************************
M=2^Nbps;  % Modulation order=1/2/4/6 for BPSK/QPSK/16QAM/64QAM
Ng=Nfft/8;
Nsym=Nfft+Ng;      % Symbol duration

% Nvc=Nfft/8;        % Nvc=0: no virtual carrier
Nvc = 0            

Nused=Nfft-Nvc;

file_name=['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '_Nframe' num2str(Nframe) '_My.dat'];
fid=fopen(file_name, 'w+');

norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM

if(meaSwitch==1)
    iBegin=0;
else
    iBegin=1;
end

for i=iBegin :length(EbN0)
   randn('state',0); rand('state',0); 
   Neb=0; Ntb=0; % Initialize the number of error/total bits
   for m=1:N_iter
       if(i==0)
           m
       end
      % Tx______________________________________________________________
      % X= randint(1,Nused*Nframe,M); % bit: integer vector
      X = randi(M,1,Nused*Nframe) -1;
      
      Xmod= qammod(X,M,0,'gray')/norms(Nbps);
      Xber=X(Nused*Nesti+1:end);        % for ber 
      
      Disp=0;
      if(Disp==1)
          figure;
          subplot(221); plot(real(X),imag(X),'*');
          subplot(223); plot(real(Xmod),imag(Xmod),'*');
          return;
      end

     % generate Nframe of OFDM symbols
      if NgType~=2
          x_GI_temp=zeros(1,Nframe*Nsym);             
          x_GI=zeros(1,Nframe*Nsym);
      elseif NgType==2
          x_GI_temp= zeros(1,Nframe*Nsym+Ng);   % Extend an OFDM symbol by Ng zeros 
          x_GI= zeros(1,Nframe*Nsym+Ng);   % Extend an OFDM symbol by Ng zeros 
      end
      kk1=[1:Nused/2]; kk2=[Nused/2+1:Nused]; kk3=1:Nfft; kk4=1:Nsym;
      for k=1:Nframe
         if Nvc~=0, X_shift= [0 Xmod(kk2) zeros(1,Nvc-1) Xmod(kk1)];
          else      X_shift= [Xmod(kk2) Xmod(kk1)];
         end
         x= ifft(X_shift);
         x_GI_temp(kk4)= guard_interval(Ng,Nfft,NgType,x);      
         
         if Nonliear==1      % ���ӷ�����PAЧӦ
            dispTemp=0;
            [x_GI(kk4),OBO] = nonlPA(x_GI_temp(kk4),dispTemp,PAbasicAmp,bq_CO,P);
         else               % �����ӷ����ԣ������л����Ŵ����ķŴ�
             x_GI(kk4)=x_GI_temp(kk4) * PAbasicAmp;             
         end
        
         kk1=kk1+Nused; kk2= kk2+Nused; kk3=kk3+Nfft; kk4=kk4+Nsym;
      end

      Disp=0;
      if(Disp==1)
          figure;
          subplot(221); plot(abs(x_GI));
          subplot(223); plot(real(x_GI),imag(x_GI),'*');
          return;
      end

      if Ch==0
        y= x_GI;  % No channel
      else  % Multipath fading channel
        channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);
        ch=zeros(1,Lch); 
        ch(Delay+1)=channel; % cir: channel impulse response
        y = conv(x_GI,ch); 
      end

      Disp=0;
      if(Disp==1)
          figure;
          subplot(321); plot(abs(x_GI));            subplot(322); plot(real(x_GI),imag(x_GI),'*');
          subplot(323); plot(abs(ch));              subplot(324); plot(real(ch),imag(ch),'*');
          subplot(325); plot(abs(y));               subplot(326); plot(real(y),imag(y),'*');
          return;
      end

      if i==0 % Only to measure the signal power for adding AWGN noise
        y1=y(1:Nframe*Nsym);
        sigPow = sigPow + y1*y1';
        continue;
      end

      % Add AWGN noise________________________________________________
      snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
      noise_mag = sqrt((10.^(-snr/10))*sigPow/2);
      y_GI = y + noise_mag*(randn(size(y))+j*randn(size(y)));

      Disp=0;
      if(Disp==1)
          figure;
          subplot(321); plot(abs(y));            subplot(322); plot(real(y),imag(y),'*');
          subplot(323); plot(abs(y_GI));         subplot(324); plot(real(y_GI),imag(y_GI),'*');
          return;
      end

      % Rx_____________________________________________________________
      kk1=(NgType==2)*Ng+[1:Nsym]; 
      kk2=1:Nfft;
      kk3=1:Nused; 
      kk4=Nused/2+Nvc+1:Nfft; 
      kk5=(Nvc~=0)+[1:Nused/2];
      
      % һframe�ڣ��� Nframe=Nesti+Ncomp��OFDM���ţ�������ͬ�������ŵ�
      % ��Ϊ�����֣�ǰ��Nesti�����ڷ������������Ncomp�����ڲ��� & BERͳ��
      % ����Nesti��OFDM���Ž��з�����ϵ���������ŵ�������Ӧ�ķ��룻
      % ����Ncomp�����Ž��л��ڷ�������Ĳ�����BERͳ��
      Ncomp=Nframe - Nesti;
      y_GI_EstiTemp=y_GI(1 : Nesti*Nsym);
      
      y_GI_Comp=y_GI(Nesti*Nsym+1 : end);
      
      dn_estiTemp = x_GI_temp(1 : Nesti*Nsym);
      dn_compTemp = x_GI_temp(Nesti*Nsym+1 : end);      
      
      % ���������  dn_estiTemp ������ PAbasicAmp����ô���Ƶõ��� ch_est Ҳ��*PAbasicAmp����֮��*
      dn_estiTemp = dn_estiTemp * PAbasicAmp;       %PA�Ļ����Ŵ���
      
      y_GI_Esti=y_GI_EstiTemp.';        %ת��
      dn_esti=dn_estiTemp.';            %ת��
      dn_comp=dn_compTemp.';
      chT=ch.';                         %ת��
      dispTemp=0;
      [bq_est,ch_est]=separa(y_GI_Esti,Nesti,dn_esti,Lbq,Lch,P,Ng,Nsym,Nfft,bq_CO,chT,dispTemp,PAbasicAmp);
      bq_est = bq_est ./ bq_est(1);     %�����һ������Ϊ���涪����һ��
      
      Hest= fft([ch_est zeros(1,Nfft-Lch)]); % Channel frequency response
      Hest_shift(kk3)= [Hest(kk4) Hest(kk5)]; 
      
      for k_Ncomp=1:Ncomp
        Y(kk2)= fft(remove_GI(Ng,Nsym,NgType,y_GI_Comp(kk1)));      %�ú���Ĳ��ֽ��в�����BERͳ��
        Y_shift=[Y(kk4) Y(kk5)];
        
        disp_0=0;
        if(disp_0==1)
            figure;
            plot(Y_shift,'.');    
            title('�������շ���');    xlabel('I');    ylabel('Q');
        end        
        
        % Ƶ�������ŵ�����
        if Ch==0
            Xmod_r(kk3) = Y_shift;
        else
            Xmod_r(kk3)= Y_shift./Hest_shift;  % Equalizer - �����ŵ�����������channel compensation
        end

        disp_0=0;
        if(disp_0==1)
            figure;
            plot(Xmod_r(kk3),'.');    
            title('Ƶ�������ŵ�������');    xlabel('I');    ylabel('Q');
        end    
        
        % ʱ��PA�ķ����Բ���
        Xtemp2=Xmod_r(kk3);
        xn2=ifft(Xtemp2);       xn2T=xn2.';
        
        dn2_0= dn_comp( (k_Ncomp-1)*Nsym+1:k_Ncomp*Nsym);
        dn2=dn2_0(Ng+1:end);
        [Dson2,GMPmode_D2] = P_onlyOdd_Dson(dn2,P);
        
        xn3=xn2T - Dson2(:,2:end) * bq_est(2:end);
        Xmod_r(kk3)=fft(xn3);
        
        disp_0=0;
        if(disp_0==1)
            figure;
            plot(Xmod_r(kk3),'.');    
            title('ʱ��PA�ķ����Բ�����');    xlabel('I');    ylabel('Q');
        end   
        
        kk1=kk1+Nsym; kk2=kk2+Nfft; kk3=kk3+Nused; kk4=kk4+Nfft; kk5=kk5+Nfft;
      end
      
      Xmod_r = Xmod_r ./ PAbasicAmp;
      
      X_r=qamdemod(Xmod_r*norms(Nbps),M,0,'gray');
      Neb=Neb+sum(sum(de2bi(X_r,Nbps)~=de2bi(Xber,Nbps)));   % ��Xber���жԱ�error����
      Ntb=Ntb+Nused*Ncomp*Nbps;  %  ��Ncomp��OFDM����
      
      if(m==N_iter)
          disp('m==N_iter');
      end
      
      if Neb>Target_neb, 
        disp('Neb>Target_neb');
        break; 
      end
      
   end      %for m=1:N_iter

   if i==0
        sigPow= sigPow/Nsym/Nframe/N_iter;
        fprintf('Signal power= %11.3e\n', sigPow);
        fprintf(fid,'%%Signal power= %11.3e\n%%EbN0[dB]       BER\n', sigPow);
   else
        Ber = Neb/Ntb;     
        fprintf('EbN0=%3d[dB], BER=%4d/%8d =%11.3e\n', EbN0(i), Neb,Ntb,Ber)
        fprintf(fid, '%d\t%11.3e\n', EbN0(i), Ber);
     
        if Ber<1e-6
            disp('Ber<1e-6');
            break;
        end
   end
end    % for i=iBegin :length(EbN0)

if (fid~=0),  fclose(fid);   end

disp('Simulation is finished');

figure
plot_ber_ch_awgn(file_name,Nbps,Ch);

return;


