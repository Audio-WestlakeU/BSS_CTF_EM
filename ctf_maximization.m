 function [A,phis,Sigmae] = ctf_maximization(x,s,Cov,Q)
% x: IxP microphone signal
% s: JxP conditional expectation of source signal
% Cov: JPxJP conditional covariance matrix of s(:)
% Q: CTF length
%
% Author: Xiaofei Li, INRIA Grenoble Rhone-Alpes
% Copyright: Perception Team, INRIA Grenoble Rhone-Alpes
%
% The algorithm is described in the papers:
% [1] Xiaofei Li, Laurent Girin and Radu Horaud. An EM Algorithm for Audio Source Separation 
% Based on the Convolutive Transfer Function. 
% IEEE Workshop on Applications of Signal Processing to Audio and Acoustics, Oct 2017
%

[J,Ps] = size(s);
[I,P] = size(x);

ACov = reshape(diag(Cov),Ps,J)';
phis = abs(s).^2+abs(ACov);

xs = zeros(I,J*Q);
ss = zeros(J*Q,J*Q);


for p=Q:Ps
    xp = x(:,p);
    sp = s(:,p:-1:p-Q+1);
    
    sp = sp.';
    sp = sp(:);
    xs = xs+xp*sp';
    %%
    
    Covs = zeros(J*Q,J*Q);
    for j1=1:J
        for j2=1:J            
            Covs((j1-1)*Q+1:j1*Q,(j2-1)*Q+1:j2*Q) = Cov((j1-1)*Ps+p:-1:(j1-1)*Ps+p-Q+1,(j2-1)*Ps+p:-1:(j2-1)*Ps+p-Q+1);            
        end
    end
    
    ss = ss+sp*sp'+Covs;
end

A = xs/ss;

Axs = A*xs';
Sigmae = (Axs+Axs'-A*ss*A')/P;

A = reshape(A,[I,Q,J]);
