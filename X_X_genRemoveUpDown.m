function [X_X] = X_X_genRemoveUpDown(X,Lch,N_start,Ldn);
% X根据Lch进行延迟与下沉，形成相对于   channel 的 卷积矩阵: X_X
% Lch：无线信道长度
Lx=length(X);
Lr=Lx+Lch-1;
X_X1 = zeros(Lr,Lch);  
for i = 1:Lch
    X_X1(i:i+Lx-1,i) = X; 
end

X_X=X_X1(N_start+1:Ldn,:);

return;

