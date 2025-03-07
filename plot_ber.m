function plot_ber(file_name,Nbps)

%MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

EbN0dB=[0:1:30];  M=2^Nbps;
ber_AWGN = ber_QAM(EbN0dB,M,'AWGN');
ber_Rayleigh = ber_QAM(EbN0dB,M,'Rayleigh');
semilogy(EbN0dB,ber_AWGN,'r:'), hold on, semilogy(EbN0dB,ber_Rayleigh,'r-')
a= load(file_name);  semilogy(a(:,1),a(:,2),'b--s');  grid on
legend('AWGN analytic','Rayleigh fading analytic', 'Simulation');
xlabel('EbN0[dB]'), ylabel('BER'); axis([a(1,1) a(end,1) 1e-8 1])