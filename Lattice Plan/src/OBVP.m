clear;clc;close all;
%% set the initial state and end state
start_pt = [0 0 0]';
start_velocity = [1 1 0]';
target_pt = [3 4 2]';
u = [0 0 0]';
T = 0;

%% Calculate T* and u*
c = [9*(norm(start_pt)^2 + norm(target_pt)^2 + 18*start_pt'*target_pt);
     -12*((target_pt' * start_pt) - (start_pt' * start_velocity));
     -3 * norm(start_velocity)^2];

%% slove OBVP
M = zeros(4, 4);
M(2, 1) = 1;
M(3, 2) = 1;
M(4, 3) = 1;
M(1, 4) = c(1);
M(2, 4) = c(2);
M(3, 4) = c(3);
e = eig(M);
for i = 1:size(e,1)
    if (real(e(i)) > 0)
        T = real(e(i));
        break;
    end
end

pos = start_pt; vel = start_velocity; position = []; velocity = [];
position = [position, pos];
velocity = [velocity, vel];
DP = target_pt - start_pt - start_velocity*T;

%% Calculate Cost
J = T + 3 * norm(DP)^2/T^3;
disp(['路径的成本是：', num2str(J)]);

%% show the path
delta_time = 0.001;
alpha = (-3*DP)/T^3;
for t = 0 : delta_time : T
% 通过积分法验证
    pos = (-1*alpha*t^3)/3 + start_velocity * t + start_pt;
    position = [position, pos];
% 通过恒加速运动的性质
%     u = -3 * DP * (t - T)/T^3;
%     pos = pos + vel*delta_time + 0.5*u*delta_time^2;
%     vel = vel + u*delta_time;
    position = [position, pos];
end

plot3(position(1,:), position(2,:), position(3,:), 'linewidth', 1.5);
xlim([0 5])
ylim([0 5])
zlim([0 5])
box on;
grid on;
axis on;
hold on;

plot3(start_pt(1), start_pt(2), start_pt(3), 'b.', 'MarkerSize', 15);
text(start_pt(1)+0.1, start_pt(2)+0.1, start_pt(3)+0.1,'start_pt','FontSize',12)
plot3(target_pt(1), target_pt(2), target_pt(3), 'r.', 'MarkerSize', 15);
text(target_pt(1)+0.1, target_pt(2)+0.1, target_pt(3)+0.1,'target_pt','FontSize',12)





