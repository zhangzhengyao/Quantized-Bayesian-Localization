close all
clear all
clc
N = 10;
L = 100;

bits =1:1:5;
ks=length(bits);
for tt=1:1:ks
% Setting the mean vector and standard deviation matrix of the target 
mu_real = [50 50]; % target mean
target_sigma = 5;% target standard deviation matrix 
for diedai=1:1:2000

target = mu_real + randn(1,2)*target_sigma; % Generating target localization

sensors =rand(N,2)*L; % Generating sensors localization 
% Generating m_i
% N : Generate N truncated Gaussian discrete random numbers
% m_low:m_i low bound
% m_up:m_i up bound
m_low= 0;
m_up= 5;
m_mean= 2;
m_simga= 10; %biao zhun cha
untruncated = makedist('Normal',m_mean,m_simga);
truncated = truncate(untruncated,m_low,m_up);
m_shengcheng = random(truncated,N,1);
%disp(m_shengcheng');

% Initialize noise variance
%sigma_n = 1; % noise variance =rand(1)*
sigma_n_snr= 0;%10.*log10(sigma_n)
sigma_n = sqrt(10.^(-(sigma_n_snr)./10));

for i = 1:N
d(1,i) = norm(sensors(i,:)-target); 

dist(1,i) = d(1,i) + m_shengcheng(i,1) + randn(1,1)*(sigma_n);

end

% Generating quantized TOA measurement---Uniform Quantization
% a= min(dist)-0.0001;
% b= max(dist)+0.0001;
% q_weishu=bits(ta) ;% quantized bits
% q_dengji=2^q_weishu;
% q_quezhi=(b-a)./(q_dengji);
% dengji=0:1:2^q_weishu-1; 
% qujiana=a:q_quezhi:a+2^q_weishu*q_quezhi; % quantized fenqu
% partition=qujiana;
% for lll=1:1:2^q_weishu
%     aaa(lll)=(qujiana(lll)+qujiana(lll+1))./2;
% end
% qq=[a,aaa,b];
% codebook=qq;
% [index,quantizeda] = quantiz(dist,partition,codebook);
%fprintf('The quantized values: %d \n',quantizeda);

training_set = dist;% sin(t*pi/500)
len = 2^bits(tt);
[partition,codebook] = lloyds(training_set,len);
[index,quantizeda] = quantiz(dist,partition,codebook);

%%% GN-ML %%%%%%
X = sensors'; % matrix for receiver positions
x = target'; % source position to be determined
sigma2 =  sigma_n.^2*ones(N,1);
r=quantizeda';%quantizeda 
iter = 50;

x=mu_real';

for i = 1:iter
    G = jacob(X, x);
    f_TOA = sqrt(sum((ones(N,1)*x'-X').^2,2));
    C_inv = diag(1./sigma2);
    x = x+inv(G'*C_inv*G)*G'*C_inv*(r-f_TOA);
    x_mlgn(:,i)=x;
    mse_ml_gn(i)= sqrt(mean((x_mlgn(:,i)'-target).^2));

end

mse_ml_gna(:,diedai)=mse_ml_gn;
mse_ml_gnb=sum(mse_ml_gna,2)./(size(mse_ml_gna,2));
end
mse_ml(tt)=mse_ml_gnb(end);
mse_ml_gnaa(tt,:)=mse_ml_gna(end,:);
end
plot(mse_ml,'-k');
save('MLxin.mat','mse_ml_gna','mse_ml_gnaa','mse_ml_gnb','mse_ml','sensors','target','dist','quantizeda');