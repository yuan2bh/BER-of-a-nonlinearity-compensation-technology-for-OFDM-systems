function [Dson,GMPmode] = P_onlyOdd_Dson(dn,P)
% GMPmode: 'GMP', 'MP', 'P'  
%
% 相应的 bq_CO=bq_CO/bq_CO(1);       %必须的；否则后面分离
%
% dn:输入矢量
% P:非线性阶数, 奇数 or 偶数
% D:输入非线性矩阵，输入dn根据P模型构建的数据矩阵；行数为dn个数，列数为bq的总数；每一行是一个dn产生
%
% 2015.11.13.
% yhl

dn_len=length(dn);
Dson_col=(P+1)/2;

Dson=zeros(dn_len,Dson_col);
for k_Dson_col=1:Dson_col
    Dson(:,k_Dson_col)=dn .* abs(dn).^(2*(k_Dson_col-1));
end
GMPmode='P';

return;

%********************************** before ********************

% if ( (G==0) && (M==0))
%     GMPmode='P';
% elseif ( (G==0) && (M~=0) )
%         GMPmode='MP';
% elseif ( (G~=0) && (M~=0) )
%             GMPmode='GMP';
% else
%     disp('不考虑 G!=0 && M==0');
%     return;
% end
% 
% % 实现公式第一项的 “输入非线性矩阵”部分
% % M==0时，实现P模型；M!=0时实现MP模型的第一项或GMP的第一项
% YO_A = zeros(len, floor((P+1)/2) * (M+1));   
% for n = 1:len
%     S=0;                % p and m 的枚举组合数
%     for p = 0:2:P-1         %实现去除偶次项
%         for m = 0:M
%             S=S+1;
%             if(n>m)     %满足条件的写入矢量，其他为0；一个n对应一个行矢量；所有n对应一个矩阵
%                         %一个p与m组合 且 满足数据存在条件的 对应一个数据项，
%                 YO_A(n,S) = dn(n-m)*abs(dn(n-m))^p;
%             end
%         end
%     end
% end
% 
% if (GMPmode=='P')
%     YO_B=[];
% else
%     % 实现公式第二项的 “输入非线性矩阵”部分
%     % G==0时，实现MP模型的第二项；G!=0时，实现GMP模型的第二项(目前其实MP模型未使用)
%     if G==0
%         YO_B = zeros(len, floor((P+1)/2-1) * (M+1) * (G+1));
%     else
%         YO_B = zeros(len, floor((P+1)/2-1) * (M+1) * G);
%     end
%     for n = 1:len
%         S=0;
%         for p = 2:2:P-1
%             for m = 0:M
%                 if G==0
%                     S=S+1;
%                     if(n>m)
%                         YO_B(n,S) = dn(n-m)*abs(dn(n-m))^p; %将输入信号写成矩阵形式
%                     end
%                 else
%                     for g = 1:G
%                         S=S+1;
%                         if(n>m && n>m+g)
%                             YO_B(n,S) = dn(n-m)*abs(dn(n-m-g))^p; %将输入信号写成矩阵形式
%                         end
%                     end
%                 end    
%             end
%         end
%     end
% end
% 
% % ouput
% %1---GMP; 2---MP;   3---P
% switch GMPmode
%     case 'GMP'
%         Dson=[YO_A,YO_B];   %矩阵高为输入数据长度；宽度为相应矢量高度,取GMP公式的第一项、第二项
%     case 'MP'
% %         D=[YO_A,YO_B];   %矩阵高为输入数据长度；宽度为相应矢量高度,取GMP公式的第一项、第二项
%         Dson=[YO_A];        %矩阵高为输入数据长度；宽度为相应矢量高度,取GMP公式的第一项
%     case 'P'
%         Dson=[YO_A];        %矩阵高为输入数据长度；宽度为相应矢量高度,只取GMP公式的第一项
%     otherwise
%         warning('GMPMode must be : GMP, MP, or P with function GMP_MP_P');
% end
% 
% 
% return;






