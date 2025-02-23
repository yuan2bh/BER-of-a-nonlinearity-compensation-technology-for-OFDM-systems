function [bq_est,ch_estT]=separaMean(y_GI_Esti,Nesti,dn_esti,Lbq,Lch,P,Ng,Nsym,Nfft,bq_CO,ch,dispTemp,PAbasicAmp)
% dn����ģ��PA������xn,���������ŵ�ch��Ϊy����AWGN������Ϊ�����ź�rn;����rn�����Ƶ��Ĺ�ʽ���Ƶõ������ŵ������������
% format long;
    
dnMem=zeros(Nfft+Ng,Nsym);
xnMem=zeros(Nfft+Ng, Nsym);
OBOVec=zeros(1,Nesti);
    
bq_estFrame=zeros(Lbq,Nesti);
ch_estFrame=zeros(Lch,Nesti);

nmse_bqFrame=zeros(Nesti,1);
mmse_chFrame=zeros(Nesti,1);

% --------------------------
bq_estCurr=zeros(Lbq,1);   
bq_estAccu=zeros(Lbq,1);
ch_estCurr=zeros(Lch,1);
ch_estAccu=zeros(Lch,1);

% --------------------------
Ldn=Nsym;
    
for kSymbol = 1:Nesti   % Each symbol of one frame add AWGN for average of simulations
    
    % ȡ��kSymbol��OFDM����
    dn= dn_esti( (kSymbol-1)*Nsym+1:kSymbol*Nsym);
    [Dson,GMPmode_D] = P_onlyOdd_Dson(dn,P);
    xn=Dson*bq_CO;     %PA output
    dnMem(:,kSymbol)=dn;
    xnMem(:,kSymbol)=xn;
    OBOVec(kSymbol)=10*log10( max(abs(xn).^2) / mean(abs(xn).^2));    
    
    % �����ź�r���������D_D����X_X�Ŀ�ʼ�����������˲̬��Ӧ������ֵ (Ng - Lch +1) / 2
    N_start = ceil( (Ng-Lch)/2 +1 );   % 256

    [D_D] = D_D_genRemUpDown(Dson,Ldn,Lbq,Lch,N_start);      % D_D ��С��Ldn - N_start  x  Lch*Lbq;

    if ((kSymbol==1)&&(dispTemp==1))
        rank_D_D=rank(D_D)  
        rank_D=rank(Dson)
    end
    
     % ȡ��kSymbol��OFDM����
     rn1= y_GI_Esti( (kSymbol-1)*Nsym+1:kSymbol*Nsym);

     rn=rn1(N_start+1:Ldn);
     chbq_estLS=pinv(D_D)*rn;     %α��, ������С���ˣ�����kroneker product : chbq

     % --------------------------
%          ȡͷһС�飬
     if kSymbol==1
         bq_estCurr =chbq_estLS(1 : Lbq ) ./ chbq_estLS(1);     %��һ��OFDM���ų��Ե�һ��Ԫ��
     else
         bq_estCurr =chbq_estLS(1 : Lbq ) ./ ch_estCurr(1);     %�ڶ���OFDM���ſ�ʼ��������һ�ε�ch���ƽ��
%              bq_estCurr=bq_estCurr ./ bq_estCurr(1);          % ΪʲôҪ���й�һ��
     end
     bq_estAccu=bq_estCurr+bq_estAccu;
%      bq_estCurr=bq_estAccu/kSymbol;

     % ��¼ÿ��SNR�µ�һ֡��OFDM���ŵ�nmseֵ
     bq_err = bq_estCurr-bq_CO;
     bq_NorSquErr = (bq_err'*bq_err) / ( (bq_CO(2:end))'*(bq_CO(2:end)) );  %ȥ����һ��ϵ��1��Ӱ�죬��������ͬ
     nmse_bqFrame(kSymbol)=10*log10(bq_NorSquErr);     % ����ƽ��ֵ

     % --------------------------
     % ����bq�Ĺ��ƣ�����PA�������Ȼ������ŵ�
     xn_est=Dson * bq_estCurr;    
     xn_xn = X_X_genRemoveUpDown(xn_est,Lch,N_start,Ldn);   % ����PA���xn�ľ������

     ch_estCurr=pinv(xn_xn) * rn;       % ���ŵ���LS����
     ch_estAccu=ch_estAccu+ch_estCurr;
%      ch_estCurr=ch_estAccu/kSymbol;

     ch_err = ch_estCurr - ch;          % ���ŵ��������
     ch_NorSquErr = (ch_err'*ch_err) / ( ch'* ch );   % ���ŵ����Ƶ����ģ�Ĺ�һ��(����ת�ã�)
     mmse_chFrame(kSymbol)=10*log10(ch_NorSquErr);

end   % end of Nesti OFDM symbols adding AWGN һ��SNR�µĲ�ͬnoiseʵ��

bq_est = mean(bq_estAccu,2);
ch_est=mean(ch_estAccu,2).';   % for���

return;

