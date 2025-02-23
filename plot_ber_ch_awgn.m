function plot_ber_ch_awgn(file_name,Nbps,Ch)

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

% Ch= 1  
% if Ch==0
%     chType='AWGN'; 
%     Target_neb=100; 
% else
%     chType='CH'; 
%     Target_neb=500*10*2*2; 
% end

if(Ch==0)                                     % AWGN下
    EbN0dB=[0:1:30];  M=2^Nbps;
    ber_AWGN = ber_QAM(EbN0dB,M,'AWGN');
%     ber_Rayleigh = ber_QAM(EbN0dB,M,'Rayleigh');
    semilogy(EbN0dB,ber_AWGN,'r:');
    hold on;  
%     hold on, semilogy(EbN0dB,ber_Rayleigh,'r-')
    a= load(file_name);  semilogy(a(:,1),a(:,2),'b--s');  grid on
    legend('AWGN analytic噪声下理论解析曲线-无非线性', '发射机存在非线性时Simulation');
    xlabel('EbN0[dB]'), ylabel('BER'); axis([a(1,1) a(end,1) 1e-6 1])
else                                             % 多径下
    EbN0dB=[0:1:30];  M=2^Nbps;
    ber_Rayleigh = ber_QAM(EbN0dB,M,'Rayleigh');
    a= load(file_name);  
    semilogy(a(:,1),a(:,2),'b--s');  grid on;  hold on; 
    semilogy(EbN0dB,ber_Rayleigh,'r-');    
    legend('Our proposed method', 'Theoretical value for Rayleigh fading channel');
    xlabel('EbN0[dB]'), ylabel('BER'); axis([a(1,1) a(end,1) 1e-4 1])
end


