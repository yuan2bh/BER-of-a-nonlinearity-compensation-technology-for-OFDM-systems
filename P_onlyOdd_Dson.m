function [Dson,GMPmode] = P_onlyOdd_Dson(dn,P)
% GMPmode: 'GMP', 'MP', 'P'  
%
% ��Ӧ�� bq_CO=bq_CO/bq_CO(1);       %����ģ�����������
%
% dn:����ʸ��
% P:�����Խ���, ���� or ż��
% D:��������Ծ�������dn����Pģ�͹��������ݾ�������Ϊdn����������Ϊbq��������ÿһ����һ��dn����
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
%     disp('������ G!=0 && M==0');
%     return;
% end
% 
% % ʵ�ֹ�ʽ��һ��� ����������Ծ��󡱲���
% % M==0ʱ��ʵ��Pģ�ͣ�M!=0ʱʵ��MPģ�͵ĵ�һ���GMP�ĵ�һ��
% YO_A = zeros(len, floor((P+1)/2) * (M+1));   
% for n = 1:len
%     S=0;                % p and m ��ö�������
%     for p = 0:2:P-1         %ʵ��ȥ��ż����
%         for m = 0:M
%             S=S+1;
%             if(n>m)     %����������д��ʸ��������Ϊ0��һ��n��Ӧһ����ʸ��������n��Ӧһ������
%                         %һ��p��m��� �� �������ݴ��������� ��Ӧһ�������
%                 YO_A(n,S) = dn(n-m)*abs(dn(n-m))^p;
%             end
%         end
%     end
% end
% 
% if (GMPmode=='P')
%     YO_B=[];
% else
%     % ʵ�ֹ�ʽ�ڶ���� ����������Ծ��󡱲���
%     % G==0ʱ��ʵ��MPģ�͵ĵڶ��G!=0ʱ��ʵ��GMPģ�͵ĵڶ���(Ŀǰ��ʵMPģ��δʹ��)
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
%                         YO_B(n,S) = dn(n-m)*abs(dn(n-m))^p; %�������ź�д�ɾ�����ʽ
%                     end
%                 else
%                     for g = 1:G
%                         S=S+1;
%                         if(n>m && n>m+g)
%                             YO_B(n,S) = dn(n-m)*abs(dn(n-m-g))^p; %�������ź�д�ɾ�����ʽ
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
%         Dson=[YO_A,YO_B];   %�����Ϊ�������ݳ��ȣ����Ϊ��Ӧʸ���߶�,ȡGMP��ʽ�ĵ�һ��ڶ���
%     case 'MP'
% %         D=[YO_A,YO_B];   %�����Ϊ�������ݳ��ȣ����Ϊ��Ӧʸ���߶�,ȡGMP��ʽ�ĵ�һ��ڶ���
%         Dson=[YO_A];        %�����Ϊ�������ݳ��ȣ����Ϊ��Ӧʸ���߶�,ȡGMP��ʽ�ĵ�һ��
%     case 'P'
%         Dson=[YO_A];        %�����Ϊ�������ݳ��ȣ����Ϊ��Ӧʸ���߶�,ֻȡGMP��ʽ�ĵ�һ��
%     otherwise
%         warning('GMPMode must be : GMP, MP, or P with function GMP_MP_P');
% end
% 
% 
% return;






