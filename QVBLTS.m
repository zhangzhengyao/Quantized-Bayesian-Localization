close all
clear all
clc
N = 10;
L = 100;

 bits=1:1:5;
 tt=length(bits);
 for ta=1:1:tt
% Setting the mean vector and standard deviation matrix of the target 
mu_real = [50 50]; % target mean---target coarse location
target_sigma = 5;% target standard deviation matrix 
for diedai=1:1:2000
target = mu_real + randn(1,2)*target_sigma; % Generating target localization

sensors =rand(N,2)*L; % Generating sensors localization 

% Generating clock offsets Delay Measurement m_i
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
%disp(m_shengcheng');

% Initialize noise variance
%sigma_n = 0.1; % noise standard variance 
sigma_n_snr= 0;%10.*log10(1./sigma_n)
sigma_n = sqrt(10.^(-(sigma_n_snr)./10));
for i = 1:N
d(1,i) = norm(sensors(i,:)-target); 
dist(1,i) = d(1,i) + m_shengcheng(i,1) + randn(1,1)*(sigma_n); % TOA measurement
 while dist(1,i) < 0
 dist(1,i) = d(1,i) + m_shengcheng(i,1) + randn(1,1)*(sigma_n); % TOA measurement
 end
end

% Generating quantized TOA measurement-----LLOYDS Optimal Quantization
training_set = dist;
len = 2^bits(ta);
[partition,codebook] = lloyds(training_set,len);
[quantized,quantizeda] = quantiz(dist,partition,codebook);


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


% Initialization parameters and hyperparameters
   
y0 = mu_real; 
 
m0 =  m_mean*ones(N,1); 
da0 =  quantizeda; 

LB = 0; % Variational lower bound Initialize variational lower bound to 0

%tol = 1e-2; % Allowable error of convergence

maxiters = 50; % Maximum number of iterations


% model solution
for it=1:maxiters

        %solutive q(x)
        sigma_y_new=1./(N./(sigma_n^2)+1./(target_sigma.^2));

        for i = 1:N


            Taylor_ex(i,:) = -(sensors(i,:)-y0(it,:))./norm(sensors(i,:)-y0(it,:));


            U1(i,:) = (sensors(i,:) + (quantizeda(1,i) - m0(i,it))*Taylor_ex(i,:));%quantizeda(1,i) m_shengcheng(i,1) m0(i,it) da0(it,i) dist(1,i), quantizeda(1,i)

        end

        U = sum(U1,1)./(sigma_n.^2) ;%./(sigma_n.^2);

        u_y_new = (U + mu_real./(target_sigma.^2)).*sigma_y_new;

        y0(it+1,:)=u_y_new;


        %solutive  q(m)
        sigma_m_new=1./(1./(sigma_n.^2)+1./(m_simga.^2)); 

        for i = 1:N

            Um1(i,:) = (quantizeda(1,i) -  norm(sensors(i,:)-y0(it+1,:))-Taylor_ex(i,:)*(y0(it+1,:)-y0(it,:))' );% quantizeda(1,i)  da0(it,i) dist(1,i) 
            Um2(i,:) =(m_mean);

        end
        
        Um =  sum(Um1,2)./(sigma_n.^2)+sum(Um2,2)./(m_simga.^2);%

        u_m_new = Um.*sigma_m_new;
        
        % Use the find function to filter the indexes of elements greater than m_up
        indice = find(u_m_new < m_low);
        
        % Use the index to get the element that satisfies the condition
        
        u_m_new(indice)=m_low;

        % Use the find function to filter out the indexes of elements greater than m_up.
        indices = find(u_m_new > m_up);
        
        % Use the index to get the element that satisfies the condition
        
        u_m_new(indices)=m_up;

       %disp(u_m_new)

       for ed=1:N
       pdf_mlow(ed) = normpdf((m_low-u_m_new(ed))./sqrt(sigma_m_new), 0, 1);% truncated_normal_a_pdf ( quantizeda(ed), 0, 1, (m_low-deee(ed))./sigma_n );
       pdf_mup(ed) = normpdf((m_up-u_m_new(ed))./sqrt(sigma_m_new), 0, 1);% truncated_normal_b_pdf ( quantizeda(ed),  deee(ed), sigma_n, m_up );
       cdf_mlow(ed) = normcdf((m_low-u_m_new(ed))./sqrt(sigma_m_new), 0, 1);% truncated_normal_a_cdf ( quantizeda(ed), deee(ed), sigma_n, m_low );
       cdf_mup(ed) = normcdf((m_up-u_m_new(ed))./sqrt(sigma_m_new), 0, 1);% truncated_normal_b_cdf ( quantizeda(ed),  deee(ed), sigma_n, m_up );
       end
       eeeddd=(pdf_mup-pdf_mlow)./(cdf_mup-cdf_mlow).*sqrt(sigma_m_new);
       u_m_newa=u_m_new-eeeddd';
       u_m_new=u_m_new-eeeddd';
        m0(:,it+1)=u_m_new;


       %% solutive  q(d)
      
        sigma_d_new=1./(1./(sigma_n^2));

        for i = 1:N

        Ud1(i,:) = ( norm(sensors(i,:)-y0(it,:)) +Taylor_ex(i,:) *(y0(it+1,:)-y0(it,:))' + m0(i,it));%m0(i,it)  m_shengcheng(i,1)

        Ud = sum(Ud1,2)./(sigma_n.^2);

        u_d_new = Ud*sigma_d_new;
        % disp(u_d_new);
        end

       for ed=1:N
       partitiona=[0,partition,100];
       condition = partitiona > quantizeda(ed);
       indices = find(condition);
       indicesa=indices(1);
       pdf_dlow(ed) = normpdf((partitiona(indicesa-1)-u_d_new(ed))./sqrt(sigma_d_new), 0, 1);% truncated_normal_a_pdf ( quantizeda(ed), 0, 1, (m_low-deee(ed))./sigma_n );
       pdf_dup(ed) = normpdf((partitiona(indicesa)-u_d_new(ed))./sqrt(sigma_d_new), 0, 1);% truncated_normal_b_pdf ( quantizeda(ed),  deee(ed), sigma_n, m_up );
       cdf_dlow(ed) = normcdf((partitiona(indicesa-1)-u_d_new(ed))./sqrt(sigma_d_new), 0, 1);% truncated_normal_a_cdf ( quantizeda(ed), deee(ed), sigma_n, m_low );
       cdf_dup(ed) = normcdf((partitiona(indicesa)-u_d_new(ed))./sqrt(sigma_d_new), 0, 1);% truncated_normal_b_cdf ( quantizeda(ed),  deee(ed), sigma_n, m_up );
       end
       eeeddda=(pdf_dup-pdf_dlow)./(cdf_dup-cdf_dlow).* sqrt(sigma_d_new);
       u_d_new=u_d_new-eeeddda';
        da0(it+1,:)=u_d_new';

        mse_y(:,it)=sqrt(mean((y0(it+1,:) -target).^2));
      
        mse_m(it,:)=sqrt(mean((m0(:,it+1) -m_shengcheng).^2));


end

mse_yaa(:,diedai,ta)=mse_y;
mse_maa(:,diedai,ta)=mse_m; % Monte Carlo results
mse_ya(:,diedai)=mse_y;
mse_ma(:,diedai)=mse_m; % Monte Carlo results

mse_yb=sum(mse_ya,2)./(size(mse_ya,2));
mse_mb=sum(mse_ma,2)./(size(mse_ma,2));

end

mse_yc(ta)=mse_yb(end);
mse_mc(ta)=mse_mb(end);

end
figure;
plot(mse_yc);
save('QVBLbits.mat','mse_yaa','mse_maa','mse_ya','mse_ma','mse_yb','mse_mb','mse_yc','mse_mc','sensors','target','dist','quantizeda');