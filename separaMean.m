function [bq_est,ch_estT]=separaMean(y_GI_Esti,Nesti,dn_esti,Lbq,Lch,P,Ng,Nsym,Nfft,bq_CO,ch,dispTemp,PAbasicAmp)
% dn经过模型PA的数据xn,经过无线信道ch后为y，与AWGN噪声后为接收信号rn;采用rn根据推导的公式估计得到无线信道与非线性因子
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
    
    % 取第kSymbol个OFDM符号
    dn= dn_esti( (kSymbol-1)*Nsym+1:kSymbol*Nsym);
    [Dson,GMPmode_D] = P_onlyOdd_Dson(dn,P);
    xn=Dson*bq_CO;     %PA output
    dnMem(:,kSymbol)=dn;
    xnMem(:,kSymbol)=xn;
    OBOVec(kSymbol)=10*log10( max(abs(xn).^2) / mean(abs(xn).^2));    
    
    % 接收信号r、卷积矩阵D_D、与X_X的开始，消除卷积的瞬态响应；经验值 (Ng - Lch +1) / 2
    N_start = ceil( (Ng-Lch)/2 +1 );   % 256

    [D_D] = D_D_genRemUpDown(Dson,Ldn,Lbq,Lch,N_start);      % D_D 大小：Ldn - N_start  x  Lch*Lbq;

    if ((kSymbol==1)&&(dispTemp==1))
        rank_D_D=rank(D_D)  
        rank_D=rank(Dson)
    end
    
     % 取第kSymbol个OFDM符号
     rn1= y_GI_Esti( (kSymbol-1)*Nsym+1:kSymbol*Nsym);

     rn=rn1(N_start+1:Ldn);
     chbq_estLS=pinv(D_D)*rn;     %伪逆, 采用最小二乘，估计kroneker product : chbq

     % --------------------------
%          取头一小组，
     if kSymbol==1
         bq_estCurr =chbq_estLS(1 : Lbq ) ./ chbq_estLS(1);     %第一个OFDM符号除以第一个元素
     else
         bq_estCurr =chbq_estLS(1 : Lbq ) ./ ch_estCurr(1);     %第二个OFDM符号开始，除以上一次的ch估计结果
%              bq_estCurr=bq_estCurr ./ bq_estCurr(1);          % 为什么要进行归一化
     end
     bq_estAccu=bq_estCurr+bq_estAccu;
%      bq_estCurr=bq_estAccu/kSymbol;

     % 记录每个SNR下的一帧内OFDM符号的nmse值
     bq_err = bq_estCurr-bq_CO;
     bq_NorSquErr = (bq_err'*bq_err) / ( (bq_CO(2:end))'*(bq_CO(2:end)) );  %去除第一个系数1的影响，与文献相同
     nmse_bqFrame(kSymbol)=10*log10(bq_NorSquErr);     % 不是平均值

     % --------------------------
     % 根据bq的估计，计算PA的输出；然后估计信道
     xn_est=Dson * bq_estCurr;    
     xn_xn = X_X_genRemoveUpDown(xn_est,Lch,N_start,Ldn);   % 产生PA输出xn的卷积矩阵

     ch_estCurr=pinv(xn_xn) * rn;       % 求信道的LS估计
     ch_estAccu=ch_estAccu+ch_estCurr;
%      ch_estCurr=ch_estAccu/kSymbol;

     ch_err = ch_estCurr - ch;          % 求信道估计误差
     ch_NorSquErr = (ch_err'*ch_err) / ( ch'* ch );   % 求信道估计的误差模的归一化(共轭转置！)
     mmse_chFrame(kSymbol)=10*log10(ch_NorSquErr);

end   % end of Nesti OFDM symbols adding AWGN 一个SNR下的不同noise实现

bq_est = mean(bq_estAccu,2);
ch_est=mean(ch_estAccu,2).';   % for输出

return;

