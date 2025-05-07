close all
clear all
clc
N = 10;
L = 100;
bits=1:1:5;
ks=length( bits);
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
m_simga= sqrt(10);%biao zhun cha
untruncated = makedist('Normal',m_mean,m_simga);
truncated = truncate(untruncated,m_low,m_up);
m_shengcheng = random(truncated,N,1);

% Initialize noise variance

sigma_n_snr=  0;%10.*log10(sigma_n)
sigma_n = sqrt(10.^(-(sigma_n_snr)./10));
for i = 1:N
d(1,i) = norm(sensors(i,:)-target); 
dist(1,i) = d(1,i) + m_shengcheng(i,1) + randn(1,1)*(sigma_n); %randn(1,1)?sigma_i,randn(1,1)*
end

% Generating quantized TOA measurement---lloyds Quantization
training_set = dist;
len = 2^bits(tt);
[partition,codebook] = lloyds(training_set,len);
[index,quantizeda] = quantiz(dist,partition,codebook);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = sensors'; 
x = target'; 
sigma2 =  sigma_n.^2*ones(N,1);
r=quantizeda';
iter = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CRLB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crb_m=N./(sigma_n.^2);
crb_ma(:,diedai)=crb_m;
crb_mc=sum(crb_ma,2)./(size(crb_ma,2));
crb_mb=inv(crb_mc);

bcrb_m=N./(sigma_n.^2)+N./(m_simga);
bcrb_ma(:,diedai)=bcrb_m;
bcrb_mc=sum(bcrb_ma,2)./(size(bcrb_ma,2));
bcrb_mb=inv(bcrb_mc);

end
crb_mccc(tt)=crb_mb;
qbcrb_mccc(tt)=bcrb_mb;

end

save('BCRBmxin.mat','qbcrb_mccc','crb_mccc','bcrb_ma','crb_ma','sensors','target','dist','quantizeda');





