function [X_X] = X_X_genRemoveUpDown(X,Lch,N_start,Ldn);
% X����Lch�����ӳ����³����γ������   channel �� �������: X_X
% Lch�������ŵ�����
Lx=length(X);
Lr=Lx+Lch-1;
X_X1 = zeros(Lr,Lch);  
for i = 1:Lch
    X_X1(i:i+Lx-1,i) = X; 
end

X_X=X_X1(N_start+1:Ldn,:);

return;

