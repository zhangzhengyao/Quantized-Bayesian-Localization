close all
clear all
clc
N = 10;
L = 100;

bits=1:1:5;
ks=length( bits);
Maxi = 2000;
target_sigma = 5;% target standard deviation matrix 
mu_real = [50 50]; % target mean
sensors =rand(N,2)*L; % Generating sensors localization 

target = mu_real + randn(1,2)*target_sigma; % Generating target localization

 for tt=1:1:ks
m_low= 0;
m_up= 5;
m_mean= 2;
m_simga= sqrt(10);%biao zhun cha
untruncated = makedist('Normal',m_mean,m_simga);
truncated = truncate(untruncated,m_low,m_up);
m_shengcheng = random(truncated,N,1);

for diedai=1:1:Maxi
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%
X = sensors'; % matrix for receiver positions
x = target'; % source position to be determined
sigma2 =  sigma_n.^2*ones(N,1);
r=quantizeda';
iter = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CRLB%%%%%%%%%%%%%%%%%%%%%%%%

crbaaa(diedai) = crlb([target' sensors']', sigma2);
crbbbb=sum(crbaaa)./Maxi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%QBCRB%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:1:N
     d(1,j) = norm(sensors(j,:)-target); 
     AZX=(mu_real-sensors(j,:))'*(mu_real-sensors(j,:))./(d(1,j));

    AZQ=abs((quantizeda(1,j)-d(1,j) -m_shengcheng(j,1)))* (1./d(1,j).*eye(2)-((mu_real-sensors(j,:))'*(mu_real-sensors(j,:)))./(d(1,j).^3));

     QQ(:,:,j) = AZX - AZQ;
     QQ1(:,:,j) = AZX;
    fai(:,j)=((mu_real-sensors(j,:))'./d(1,j))*1./(sigma_n.^2);
    xigema=1./(m_simga.^2)*eye(N);
    sanjiao=1./(sigma_n.^2)*eye(N);

 end

 Omaxtri = sum(QQ,3);
 Omaxtri1 = sum(QQ1,3);
 qbcrbaa=((Omaxtri*(1./(sigma_n.^2)))+(1./target_sigma.^2)*eye(2));
 qbcrbaa1=((Omaxtri1*(1./(sigma_n.^2)))+(1./target_sigma.^2)*eye(2));
AAAA=inv(qbcrbaa-fai*inv(sanjiao +xigema)*fai');
AAAA1=inv(qbcrbaa1-fai*inv(sanjiao +xigema)*fai');
AAAB=trace(AAAA);
AAAB1=trace(AAAA1);
 qbcrba(:,diedai)=AAAB;
 qbcrba1(:,diedai)=AAAB1;

end
qbcrbbQ=sum(qbcrba,2)./Maxi;
qbcrbee(tt)=qbcrbbQ;
%qbcrbbQ1=sum(qbcrba1,2)./Maxi;
%qbcrbee1(tt)=qbcrbbQ1;
crbccc(tt)=crbbbb;
end

plot(qbcrbee);
hold on ;
save('BCRBxin.mat','qbcrba','qbcrbee','crbccc','sensors','target','dist','quantizeda');




