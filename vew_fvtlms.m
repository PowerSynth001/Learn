function [y, e, L_hist, sigma_hist] = vew_fvtlms(d, x, max_L, L_init)
% 输入参数:
% d: 期望信号
% x: 输入信号
% max_L: 最大允许阶数
% L_init: 初始阶数
%
% 输出参数:
% y: 滤波器输出
% e: 误差信号
% L_hist: 阶数变化历史
% sigma_hist: 误差宽度变化历史

    N = length(x);
    w = zeros(max_L, 1);    % 权系数初始化
    L_frac = L_init;        % 分数阶数初始化
    
    y = zeros(N, 1);
    e = zeros(N, 1);
    L_hist = zeros(N, 1);
    sigma_hist = zeros(N, 1);
    
    Delta = 2;              % 片段长度差，用于计算片段误差 
    rho = 0.99;             % 平滑参数 (接近于1) 
    gamma = 0.01;           % 调节常数 
    sigma = 0.1;            % 初始误差宽度 (需设为非零小值)
    
    mu_w = 0.005;           % 权值更新步长
    mu_L = 0.05;            % 阶数更新步长
    leak = 0.01;            % 阶数泄漏因子 (防止漂移)

    for n = max_L:N
        % 1. 获取当前整数阶数
        L_curr = round(L_frac);
        % 限制阶数在合法范围内
        L_curr = max(Delta + 2, min(max_L, L_curr));
        
        % 2. 计算完整误差 (Full Error)
        x_full = x(n:-1:n-L_curr+1);
        y_full = w(1:L_curr)' * x_full;
        e_full = d(n) - y_full;
        
        % 保存输出
        y(n) = y_full;
        e(n) = e_full;
        L_hist(n) = L_frac;
        sigma_hist(n) = sigma;
        
        % 3. 计算片段误差 (Segment Error)
        % 利用 L(n) - Delta 长度的滤波器输出 
        L_seg = L_curr - Delta;
        x_seg = x(n:-1:n-L_seg+1);
        y_seg = w(1:L_seg)' * x_seg;
        e_seg = d(n) - y_seg;
        
        % 4. 计算误差平方差 BE(n)
        BE = abs(e_full^2 - e_seg^2);
        
        % 5. 动态更新误差宽度 sigma(n)
        sigma = rho * sigma + gamma * BE;
        
        % 6. 更新分数阶数 L(n)
        % 使用动态的 sigma 作为分母或阈值来更新阶数
        % 这里采用标准分数阶数更新形式，将固定误差宽度替换为动态 sigma
        % 逻辑：当 e^2 > sigma 时阶数增加，反之减小
        
        gradient_L = leak * L_frac - (e_full^2) / (sigma + eps); 
        L_frac = L_frac - mu_L * gradient_L;
        
        % 阶数限幅
        L_frac = max(Delta + 2, min(max_L, L_frac));
        
        % 7. 更新权系数 (LMS)
        w(1:L_curr) = w(1:L_curr) + mu_w * e_full * x_full;
        if L_curr < max_L
            w(L_curr+1:end) = 0; % 清除无效阶数的权值
        end
    end
end