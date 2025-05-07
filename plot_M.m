
load('BCRBmxin.mat');
load("QVBLbits.mat");
figure;

plot(mse_mc,'r--o', 'LineWidth', 3);
 hold on;
 grid on;
 plot(crb_mccc, 'b-x', 'LineWidth', 3);
 plot(qbcrb_mccc,"m-.", 'LineWidth', 3); 
 xlabel('Quantized Bits');
 ylabel('RMSE of Clock Offset');

set(gca,'XTick',1:1:5);
 for ta=1:1:5
y = (mse_maa(50,:,ta));
 AX=size(y,2);
error(ta,:) = std(y');
yy=(mse_mc)';
end
errorbar(yy,error ,'.','Color', 'red','LineWidth', 2);
ylim([-2 6]);
set(gca,'YTick',-2:2:6); 
 legend('QVBL-TS','CRB','QBCRB');%legend('VBI','LS','CRB','QCRB');,'LLS','BCRB'