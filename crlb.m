function v = crlb(X, sigma2)
% CRLB
% --------------------------------
% v = crlb(X, sigma2)
% v = sum of the CRLB diagonal elements
% X = first row is the source position; second to last rows are receiver
% positions
% sigma2 = noise variance vector
% 
k = size(X, 1);
M = k - 1;
d = sqrt(sum((X(2:end, :) - ones(M, 1)*X(1, :)).^2,2));
N = size(X,2);
D = (ones(M, 1)*X(1, :) - X(2:end, :))./(d*ones(1, N));
FIM = D'*diag(1./sigma2)*D;
H = inv(FIM);
v = trace(H);
