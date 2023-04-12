clc;clear;close all;
path = ginput() * 100.0;
% path = [8.74769797421731,60.8644859813084;
%         22.5598526703499,30.9579439252336;
%         36.5561694290976,67.8738317757009;
%         49.4475138121547,34.9299065420561;
%         69.1528545119705,65.0700934579439;
%         83.5174953959484,34.2289719626168];

n_order       = 7;% order of poly
n_seg         = size(path,1)-1;% segment number
n_poly_perseg = (n_order+1); % coef number of perseg

ts = zeros(n_seg, 1);
% calculate time distribution in proportion to distance between 2 points
dist     = zeros(n_seg, 1);
dist_sum = 0;
T        = 25;
t_sum    = 0;

for i = 1:n_seg
    dist(i) = sqrt((path(i+1, 1)-path(i, 1))^2 + (path(i+1, 2) - path(i, 2))^2);
    dist_sum = dist_sum+dist(i);
end
for i = 1:n_seg-1
    ts(i) = dist(i)/dist_sum*T;
    t_sum = t_sum+ts(i);
end
ts(n_seg) = T - t_sum;

% or you can simply set all time distribution as 1
% for i = 1:n_seg
%     ts(i) = 1.0;
% end

poly_coef_x = MinimumSnapQPSolver(path(:, 1), ts, n_seg, n_order);
poly_coef_y = MinimumSnapQPSolver(path(:, 2), ts, n_seg, n_order);

% display the trajectory
X_n = [];
Y_n = [];
k = 1;
tstep = 0.01;
for i=0:n_seg-1
    %#####################################################
    % STEP 3: get the coefficients of i-th segment of both x-axis
    % and y-axis
    for t = 0:tstep:ts(i+1)
        Pxi = flipud(poly_coef_x(i*(n_order+1)+1:(i+1)*(n_order+1),1));
        Pyi = flipud(poly_coef_y(i*(n_order+1)+1:(i+1)*(n_order+1),1));
        X_n(k)  = polyval(Pxi, t);
        Y_n(k)  = polyval(Pyi, t);
        k = k + 1;
    end
end
 
plot(X_n, Y_n , 'g', 'LineWidth', 2);
hold on;
grid on;
scatter(path(1:size(path, 1), 1), path(1:size(path, 1), 2));
plot(path(:,1), path(:,2), 'r--');

function poly_coef = MinimumSnapQPSolver(waypoints, ts, n_seg, n_order)
    start_cond = [waypoints(1), 0, 0, 0];   % p, v, a, j
    end_cond   = [waypoints(end), 0, 0, 0]; % p, v, a, j
    %#####################################################
    % STEP 1: compute Q of p'Qp
    Q = getQ(n_seg, n_order, ts, 4);
    %#####################################################
    % STEP 2: compute Aeq and beq
    [Aeq, beq] = getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond);
    f = zeros(size(Q,1),1);
    poly_coef = quadprog(Q,f,[],[],Aeq, beq);
end