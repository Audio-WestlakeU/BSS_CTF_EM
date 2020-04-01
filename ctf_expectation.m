function [s,Cov] = ctf_expectation(x,A,phis,Sigmae)
%
% x: IxP microphone signals
% A: IxQxJ CTF coefficients
% phis: JxPs variance of source signal
% Sigmae: IxI covariance of noise signal
%
% Author: Xiaofei Li, INRIA Grenoble Rhone-Alpes
% Copyright: Perception Team, INRIA Grenoble Rhone-Alpes
%
% The algorithm is described in the papers:
% [1] Xiaofei Li, Laurent Girin and Radu Horaud. An EM Algorithm for Audio Source Separation 
% Based on the Convolutive Transfer Function. 
% IEEE Workshop on Applications of Signal Processing to Audio and Acoustics, Oct 2017
%

[I,P] = size(x);
[~,Q,J] = size(A);
Ps = P-Q+1;

%%
Ac = zeros(I*P,J*Ps); 
for i=1:I
    for j = 1:J        
        Aij = squeeze(A(i,:,j));
        for p = 1:Ps
            endp = min(P,p+Q-1);
            Ac((i-1)*P+p:(i-1)*P+endp,(j-1)*Ps+p) = Aij(1:endp-p+1);
        end
    end
end

phis = phis';
Psisinv = diag(1./phis(:));

Psie = zeros(I*P,I*P);
for i1 = 1:I
    for i2 = 1:I
        Psie((i1-1)*P+1:i1*P,(i2-1)*P+1:i2*P) = Sigmae(i1,i2)*eye(P);
    end
end

x = x.';
x = x(:);

Cov = inv(Ac'*(Psie\Ac)+Psisinv);

shat = Cov*(Ac'*(Psie\x));
s = reshape(shat,[Ps,J]).';



















