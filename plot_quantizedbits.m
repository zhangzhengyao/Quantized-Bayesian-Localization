clc;
clear all;
bits=1:1:5;
load("LSxin.mat");
load("MLxin.mat");
load("QVBLbits.mat");
load("YANxin.mat");
load('BCRBxin.mat');

figure;
plot(bits,mse_yc,'r--o', 'LineWidth', 3);
 hold on;
 grid on;
 plot(bits, mse_ls, 's-g', 'LineWidth', 3);
 plot(bits, mse_ml, '^--b', 'LineWidth', 3);
 plot(bits, mse_ycyan, 'x-c', 'LineWidth', 3);
 plot(bits,crbccc,"b-x", 'LineWidth', 3);
 plot(bits,qbcrbee,"m-.", 'LineWidth', 3); 
 xlabel('Quantized Bits');
 ylabel('RMSE of Target Localization');
set(gca,'XTick',1:1:5);

for ta=1:1:5
y = mse_yaa(50,:,ta);
 AX=size(y,2);
error(ta,:) = std(y');
yy=(mse_yc)';
end
errorbar(yy,error ,'.', 'Color', 'red','LineWidth', 2);
hold on;
for ta=1:1:5
y = mse_ls_nraa(ta,:);
 AX=size(y,2);
error(ta,:) = std(y');
yy=(mse_ls)';
end
errorbar(yy,error ,'.','Color', 'green','LineWidth', 2);

hold on;
for ta=1:1:5
y = mse_ml_gnaa(ta,:);
 AX=size(y,2);
error(ta,:) = std(y');
yy=(mse_ml)';
end
errorbar(yy,error ,'.','Color', 'blue','LineWidth', 2);

hold on;
for ta=1:1:5
y = mse__ycyana(ta,:);
 AX=size(y,2);
error(ta,:) = std(y');
yy=(mse_ycyan)';
end
errorbar(yy,error ,'.','Color', '[0,1,1]','LineWidth', 2);

% AX=size(y,1);
%  XTick=1:1:AX;
% % set(gca, 'XTick',XTick);
% set(gca, 'XTick',XTick);
% xlabel('Iter Number');
% ylabel('Standard Deviations Error Bar');
ylim([-0.5 10]);
set(gca,'YTick',-1:1:10); 
legend('QVBL-TS','LLS','GN-ML','SDP-Yan','CRB','QBCRB'); 
