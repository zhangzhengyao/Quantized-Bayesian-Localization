close all
clear all
clc
N = 10;
L = 100;
%load("b.mat");
bits=1:1:5;
ks=length( bits);
for tt=1:1:ks
% Setting the mean vector and standard deviation matrix of the target 
mu_real = [50 50]; % target mean
target_sigma = 5;% target standard deviation matrix 


target = mu_real + randn(1,2)*target_sigma; % Generating target localization



sensors =rand(N,2)*L; % Generating sensors localization 
for diedai=1:1:2000
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
%T_u=1;
% Initialize noise variance
%sigma_n = 1; % noise variance =rand(1)*
sigma_n_snr=  0;%10.*log10(sigma_n)
sigma_n = sqrt(10.^(-(sigma_n_snr)./10));
% dB = 30; % SNR (in dB) which is defined as the mean of squared distance over noise variance
for i = 1:N
d(1,i) = norm(sensors(i,:)-target); 
% sigma_n(1,i) = sqrt(d(1,i).^2/10^(dB/10)); % sigma2--square of sigma, here we use: SNR_dB = 10log(d^2/sigma^2) 0.01*ones(L,1);%
dist(1,i) = d(1,i) + m_shengcheng(i,1) + randn(1,1)*(sigma_n); %randn(1,1)?sigma_i,randn(1,1)*
% Taylor_ex(i,:) = -(sensors(i,:)-target)/d(1,i);
end
%m_shengcheng = unifrnd(0,T_u, [N, 1]); % Generating clock deviations
%Generating measured values
%Generating quantized values
% a=max(dist)+0.001;
% b=min(dist)-0.001;
% q_weishu=3;%quantized bits
% q_dengji=2^q_weishu;
% q_quezhi=(a-b)./(q_dengji);
% dengji=0:1:2^q_weishu+1;
% partition=[b,b+q_quezhi,b+2*q_quezhi,b+3*q_quezhi,b+4*q_quezhi,b+5*q_quezhi,b+6*q_quezhi,b+7*q_quezhi,b+8*q_quezhi];
% des=[0,partition,a+1];
% for h=1:1:N
%     desa(h)=(des(h)+des(h+1))./2;
% end
% codebook=desa;
% [index,quantizeda] = quantiz(dist,partition,codebook);
% training_set = dist;% sin(t*pi/500)
% 
% len = 2^bits(tt);
% [partition,codebook] = lloyds(training_set,len);
% % signal = cos(t*pi/250);
% [quantized,quantizeda] = quantiz(dist,partition,codebook);

a=min(dist)-0.00001; %quantizaed min value
b=max(dist)+0.00001; %quantizaed max value

q_dengji=2^bits(tt);
q_quezhi=(b-a)./(q_dengji);
dengji=0:1:2^bits(tt)-1; 
qujiana=a+q_quezhi:q_quezhi:a+(2^bits(tt)-1)*q_quezhi; % quantized fenqu
partition=qujiana;
qujianb=[a,qujiana,b];
for lll=1:1:(2^bits(tt))
    aaaaaa(lll)=(qujianb(lll)+qujianb(lll+1))./2;
end
%aaaaaa
codebook=aaaaaa;
[index,quantizeda] = quantiz(dist,partition,codebook);


% % fprintf('The quantized values: %d \n',quantized)
% 
% qq=[(0+b)/2,(b+b+q_quezhi)/2,(b+q_quezhi+b+2*q_quezhi)/2,(b+2*q_quezhi+b+3*q_quezhi)/2,(b+3*q_quezhi+b+4*q_quezhi)/2,(b+4*q_quezhi+b+5*q_quezhi)/2,(b+5*q_quezhi+b+6*q_quezhi)/2,(b+6*q_quezhi+b+7*q_quezhi)/2,(b+7*q_quezhi+b+8*q_quezhi)/2,2*
% partition=[b,b+q_quezhi,b+2*q_quezhi,b+3*q_quezhi,b+4*q_quezhi,b+5*q_quezhi,b+6*q_quezhi,b+7*q_quezhi,b+8*q_quezhi];
% codebook=qq;
% [indexa,quantizeda] = quantiz(dist,partition,codebook);
% % fprintf('The quantized values: %d \n',quantizeda);
% Generating quantized values
%t = dist;


% % training_set = dist;% sin(t*pi/500)
% % len = 8;
% % [partition,codebook] = lloyds(training_set,len);
% % % signal = cos(t*pi/250);
% % [index,quantizeda] = quantiz(dist,partition,codebook);
% % plot(dist,quants,'.');
% 
% % Generating quantized values
% a=min(dist)-0.0001;
% b=max(dist)+0.0001;
% q_weishu=3 ;%quantized bits
% % ka=length(q_weishu);
% % for k=1:1:ka
% q_dengji=2^q_weishu(k);
% q_quezhi=(b-a)./(q_dengji);
% dengji=0:1:2^q_weishu(k)-1; 
% qujiana=a:q_quezhi:a+2^q_weishu(k)*q_quezhi; % quantized fenqu
% partition=qujiana;
% for lll=1:1:2^q_weishu(k)
%     aaa(lll)=(qujiana(lll)+qujiana(lll+1))./2;
%     % cc=1:1:2^q_weishu(k);
%     % bb(cc)=b+cc*q_quezhi;
% end
% qq=[a,aaa,b];
% codebook=qq;
% [index,quantizeda] = quantiz(dist,partition,codebook);
% % end
% % %fprintf('The quantized values: %d \n',quantizeda);


% training_set = dist;% sin(t*pi/500)
% len = 2^bits(tt);
% [partition,codebook] = lloyds(training_set,len);
% % signal = cos(t*pi/250);
% [index,quantizeda] = quantiz(dist,partition,codebook);

%%% NR-LS + GN-ML %%%%%%
X = sensors'; % matrix for receiver positions
x = target'; % source position to be determined
% x = [2,3]'; % source position to be determined
%L = size(X,2); % number of receivers
% d = (sqrt(sum((x*ones(1,L)-X).^2,1))).'; %noise-free ranges 
%dB = 30; % SNR (in dB) which is defined as the mean of squared distance over noise variance
sigma2 =  sigma_n.^2*ones(N,1);%d.^2/10^(dB/10); % sigma2--square of sigma, here we use: SNR_dB = 10log(d^2/sigma^2)    0.01*ones(L,1);%
% SNR_dB = 10*log10(0.01);
%r = d + m_shengcheng + sqrt(sigma2);%randn(L,1).*   m_shengcheng + 
r=quantizeda';%quantizeda??dist
iter = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CRLB%%%%%%%%%%%%%%%%%%%%%%%%
 % for j=1:1:N
 %     aaa(j)=(ones(1,2)-sensors(j,:))*(ones(1,2)-sensors(j,:))'./((sensors(j,:)-target)*(sensors(j,:)-target)'*sigma_n.^2);
 %     bbb(j)=(ones(1,2)-sensors(j,:))*(ones(1,2)-sensors(j,:))'*(dist(1,j) - d(1,j) - m_shengcheng(j,1))./((sensors(j,:)-target)*(sensors(j,:)-target)'*sqrt((sensors(j,:)-target)*(sensors(j,:)-target)')*sigma_n.^2);
 % end
for j=1:1:N
     aaa(j)=-(target-sensors(j,:))*(target-sensors(j,:))'./((sensors(j,:)-target)*(sensors(j,:)-target)');%*sigma_n.^2
     %bbbb(j)=(norm(sensors(j,:)-target)-(target-sensors(j,:))*(target-sensors(j,:))'./(norm(sensors(j,:)-target)))./((target-sensors(j,:))*(target-sensors(j,:))')*(quantizeda(1,j) - d(1,j) - m_shengcheng(j,1));
 end
 crb=-1./(sum(aaa)./(sigma_n.^2));
 crba(:,diedai)=crb;
 crbb=sum(crba,2)./(size(crba,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%QCRB%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:1:N
     aaaa(j)=-(target-sensors(j,:))*(target-sensors(j,:))'./((sensors(j,:)-target)*(sensors(j,:)-target)');%*sigma_n.^2
     bbbb(j)=(norm(sensors(j,:)-target)-(target-sensors(j,:))*(target-sensors(j,:))'./(norm(sensors(j,:)-target)))./((target-sensors(j,:))*(target-sensors(j,:))')*(quantizeda(1,j) - d(1,j) - m_shengcheng(j,1));
 end
 qcrb=-1./((sum(aaaa)+sum(bbbb))./(sigma_n.^2));
 qcrba(:,diedai)=qcrb;
 qcrbb=sum(qcrba,2)./(size(qcrba,2));

 

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%QBCRB%%%%%%%%%%%%%%%%%%%%%%%%%
%  for j=1:1:N
%      qaaaa(j)=(mu_real-sensors(j,:))*(mu_real-sensors(j,:))'./((sensors(j,:)-mu_real)*(sensors(j,:)-mu_real)');%*sigma_n.^2
%      qbbbb(j)=-(quantizeda(1,j)-norm(sensors(j,:)-mu_real)-m_mean)* (norm(sensors(j,:)-mu_real)-(mu_real-sensors(j,:))*(mu_real-sensors(j,:))'./(norm(sensors(j,:)-mu_real)))./((mu_real-sensors(j,:))*(mu_real-sensors(j,:))');
%      %qbbbb(j)=-(qiuantizeda(1,j)-norm(sensors(j,:)-mu_real)-m_mean);
%      %(norm(sensors(j,:)-mu_real)-(mu_real-sensors(j,:))*(mu_real-sensors(j,:))'./(norm(sensors(j,:)-mu_real)))./((mu_real-sensors(j,:))*(mu_real-sensors(j,:))')*(quantizeda(1,j) - norm(sensors(j,:)-mu_real) - m_mean);
%  end
%  qbcrb=1./((sum(qaaaa)+sum(qbbbb))./(sigma_n.^2)+target_sigma);
%  qbcrba(:,diedai)=qbcrb;
%  qbcrbb=sum(qbcrba,2)./(size(qbcrba,2));

% mse_ml_gna(:,diedai)=mse_ml_gn;
% mse_ml_gnb=sum(mse_ml_gna,2)./(size(mse_ml_gna,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%QBCRB%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:1:N
     d(1,j) = norm(sensors(j,:)-mu_real); 
     AZX=(mu_real-sensors(j,:))'*(mu_real-sensors(j,:))./(d(1,j));
     qqw(:,j)=AZX(:,1);%*sigma_n.^2
     qqwa(:,j)=AZX(:,2);%*sigma_n.^2
     %qaaaa=
    AZQ=-(quantizeda(1,j)-norm(sensors(j,:)-mu_real)-m_mean)* (1./d(1,j).*eye(2)-((mu_real-sensors(j,:))'*(mu_real-sensors(j,:)))./(d(1,j).^3));
     qqe(:,j)=AZQ(:,1);%*sigma_n.^2
     qqea(:,j)=AZQ(:,2);%*sigma_n.^2
     %qbbbb(j,:)=
    fai(:,j)=((mu_real-sensors(j,:))'./d(1,j))./(sigma_n.^2);
    xigema=1./m_simga.^2*eye(N);
    sanjiao=1./sigma_n.^2*eye(N);
%qbbbb(j,:)=
     %qbbbb(j)=-(qiuantizeda(1,j)-norm(sensors(j,:)-mu_real)-m_mean);
     %(norm(sensors(j,:)-mu_real)-(mu_real-sensors(j,:))*(mu_real-sensors(j,:))'./(norm(sensors(j,:)-mu_real)))./((mu_real-sensors(j,:))*(mu_real-sensors(j,:))')*(quantizeda(1,j) - norm(sensors(j,:)-mu_real) - m_mean);
 end
 as=sum(qqw,2);
 asa=sum(qqwa,2);
 asae=[as,asa];
  ae=sum(qqe,2);
 aea=sum(qqea,2);
 asee=[ae,aea];
 qbcrbaa=(((asee+asae)./(sigma_n.^2))+(1./target_sigma.^2)*eye(2));
AAAA=inv(qbcrbaa-fai*inv(xigema+sanjiao)*fai');
AAAB=trace(AAAA);
 qbcrba(:,diedai)=AAAB;
qbcrbb=sum(qbcrba,2)./(size(qbcrba,2));
end
%crbccc(tt)=crbbbb;
crbddd(tt)=crbb;
qcrbccc(tt)=qcrbb;
qbcrbee(tt)=sqrt(qbcrbb);

end
% save('mycrbdata_bits.mat','crbccc','qcrbccc','qbcrbee');
% plot(mse_ml_gnb,'-k');
%save('mymldata_bits.mat', 'mse_ml', '-append');
plot(qbcrbee);
%hold on ;
%plot(crbccc);





