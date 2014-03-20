% Function MRCPS 
% MRCPS is consensus clustering method for stratifying cancer 
% patients, based on both molecular and clinical data. 
% Input: S - Categorical coassociation matrix (N-by-N); N-Number of samples
%        W - Numerical molecular density matrix (N-by-N);
%        k - Number of sample clusters;
%        MaxIter - Maximum number of interations;
%        epsilon - Allowed error
%        U0 - Initial U         
% Output:U - Binary membership matrix (N-by-k);
%        D - Diagonal matrix indicating the sizes of the clusters (k-by-k);
%        error - The Frobium norm error.


function [U D error]=MorCPS(S,W,k,MaxIter,epsilon,U0)

% Initialization
t=1;
delta=realmax;
N=size(W,1);
% U0=rand(N,k);
D0=rand(k,k);
pre_delta=0;
while((t<MaxIter)&&(delta>epsilon))%&&(abs(delta-pre_delta)>=10e-7)
    t=t+1;
    pre_delta=delta;
    % update U
    temp1=(W.*S)*U0;
    temp2=(W.*(U0*D0*U0'))*U0;
    for i=1:N
        for j=1:k
            if( (temp1(i,j)-0>=10e-20) && (temp2(i,j)-0>=10e-20))
                U0(i,j)=U0(i,j)*sqrt(temp1(i,j)/temp2(i,j));
            end
        end
    end
    % update D
    temp1=U0'*temp1;
    temp2=U0'*temp2;
    for i=1:k
        for j=1:k
            D0(i,j)=D0(i,j)*sqrt(temp1(i,j)/temp2(i,j));
        end
    end
    % Update delta
    temp3=(W.*S)-(W.*(U0*D0*U0'));
    delta=norm(temp3,'fro');
    t,
    delta,
end
error=delta;
U=U0;
D=D0;
