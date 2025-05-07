load("y.mat");
figure;  

% N = 10;
% L = 100;
% % Setting the mean vector and standard deviation matrix of the target 
% mu_real = [50 50]; % target mean---target coarse location
% target_sigma = 5;% target standard deviation matrix 
% 
% target = mu_real + randn(1,2)*target_sigma; % Generating target localization
% 
% %sensors =rand(N,2)*L; % Generating sensors localization 
% sensors=[5,5;5,90;90,5;90,91;25,80;5,50;90,50;75,90;50,10;50,80];

hold on;  
markerSize = 10;
for i=1:1:10
h1=plot(sensors(i, 1),sensors(i, 2),'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 'r');  
% plot(sensors(2, 1),sensors(2, 2),'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 'r'); 
% plot(sensors(3, 1),sensors(3, 2),'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 'r'); 
% plot(sensors(4, 1),sensors(4, 2),'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 'r'); 
% plot(sensors(5, 1),sensors(5, 2),'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 'r'); 
% legend('Sensor Node i');
end
hold off;

hold on;
h2=plot(target(1),target(2),  'O', 'MarkerSize', markerSize + 2, 'MarkerFaceColor', 'b');  
 
% title('???????????????????');  
xlabel('X/m');  
ylabel('Y/m');  
legend([h1, h2], 'Sensor Node', 'Target Node');
%legend('Sensor Node 1','Sensor Node 2','Sensor Node 3','Sensor Node 4','Sensor Node 5','Sensor Node 6','Sensor Node 7','Sensor Node 8','Sensor Node 9','Sensor Node 10','Target Node');

grid on;  