
close all
clear all
clc

N = 10;
L = 100;

bits =1:1:5;
ks=length(bits);
for tt=1:1:ks
% Setting the mean vector and standard deviation matrix of the target 
mu_real = [50 50]'; % target mean
 target_sigma = 5;% target standard deviation matrix 

for diedai=1:1:2000
target = mu_real + randn(2,1)*target_sigma; % Generating target localization
sensors =rand(2,N)*L; % Generating sensors localization 

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
%disp(m_shengcheng');

% Initialize noise variance
sigma_n_snr= 0;%10.*log10(sigma_n)
sigma_n = sqrt(10.^(-(sigma_n_snr)./10));

for i = 1:N
d(i,1) = norm(sensors(:,i)-target); 
dist(i,1) = d(i,1) + m_shengcheng(i,1) + randn(1,1)*(sigma_n); 
end

training_set = dist;

len = 2^bits(tt);

[partition,codebook] = lloyds(training_set,len);

[quantized,quantizeda] = quantiz(dist,partition,codebook);

quantized_value= dec2bin(quantized);

m_shengcheng = 0*ones(N,1);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',500);

%fun = @(x)(quantizeda'-m_shengcheng-sqrt((sensors-x).^2+(sensors-x).^2))'*((quantizeda'-m_shengcheng-sqrt((sensors-x).^2+(sensors-x).^2)));

x0 = [50,50]';
[x,fval,exitflag,output] = fminsearch(@(x)(objectivefcn11(N,quantizeda,m_shengcheng,sensors,x)),x0,options) ;

mse_y=sqrt(mean((x -target).^2));
mse_ya(:,diedai)=mse_y;
mse_yb=sum(mse_ya,2)./(size(mse_ya,2));

end  

mse_ycyan(tt)=mse_yb;
mse__ycyana(tt,:)=mse_ya(end,:);

end

figure;
plot(mse_ycyan);

save('YANxin.mat','mse_ya','mse_yb','mse_ycyan','mse__ycyana','sensors','target','dist','quantizeda');