function [D_D] = D_D_genRemUpDown(D,Ldn,Lbq,Lch,N_start);
% D����Lch�����ӳ����³����γ������   channel and ���������ӵ�kronech���� �������: D_D
% rn: �����ź�ʸ��
% D: ����ÿһ����ʸ���� �������ţ�1�������� ����GMP������ģ�ͽṹ����
%      һ�����շ��Ŷ�ӦD��һ��
% Ldn: ����������ų���
% Lbq: PA������ϵ���ܸ���
% Lch�������ŵ�����
Lr=Ldn+Lch-1;
% Lr=length(rn);
D_D1 = zeros(Lr,Lbq*Lch);  
for i = 1:Lch
    D_D1(i:i+Ldn-1,(i-1)*Lbq+1:i*Lbq) = D; 
end

D_D=D_D1(N_start+1:Ldn, :);   % D_D ��С��Ldn - N_start  x  Lch*Lbq;




