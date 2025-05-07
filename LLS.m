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
m_simga= 10;%biao zhun cha
untruncated = makedist('Normal',m_mean,m_simga);
truncated = truncate(untruncated,m_low,m_up);
m_shengcheng = random(truncated,N,1);

% Initialize noise variance
%sigma_n = 1; % noise variance =rand(1)*
sigma_n_snr= 0;%10.*log10(sigma_n)
sigma_n = sqrt(10.^(-(sigma_n_snr)./10));
for i = 1:N
d(1,i) = norm(sensors(i,:)-target); 
dist(1,i) = d(1,i) + m_shengcheng(i,1) + randn(1,1)*(sigma_n); 
end

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

%%% LLS%%%%%%
sensors_x=sensors(:,1);
sensors_y=sensors(:,2);
for i=1:N-1
A(i,:)=[-2*(sensors_x(i+1)-sensors_x(1)),-2*(sensors_y(i+1)-sensors_y(1))];
b(i,:)=quantizeda(i+1).^2-sensors_x(i+1).^2-sensors_y(i+1).^2-quantizeda(1).^2+sensors_x(1).^2+sensors_y(1).^2;

end
Theta=inv(A.'*A)*A.'*b;
X=Theta(1);
Y=Theta(2);

mse_ls_nr= sqrt(mean((Theta'-target).^2));

mse_ls_nra(:,diedai)=mse_ls_nr;
mse_ls_nrb=sum(mse_ls_nra,2)./(size(mse_ls_nra,2));
end
mse_ls(tt)=mse_ls_nrb(end);
mse_ls_nraa(tt,:)=mse_ls_nra;
end
plot(mse_ls,'-b');
save('LSxin.mat','mse_ls_nr','mse_ls_nra','mse_ls_nraa','mse_ls','sensors','target','dist','quantizeda');