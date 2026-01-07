function [y, e, w] = htvss_lms(d, x, L)
% 输入参数:
% d: 期望信号 (列向量)
% x: 输入信号 (列向量)
% L: 滤波器阶数
%
% 输出参数:
% y: 滤波器输出
% e: 误差信号
% w: 最终滤波器权系数

    N = length(x); % 信号长度
    w = zeros(L, 1); % 初始化权系数
    y = zeros(N, 1); % 初始化输出信号
    e = zeros(N, 1); % 初始化误差信号
    
    % 算法相关变量初始化
    beta = 0;   % 误差信号累积值的平均值
    p = 0;      % 误差信号累积值
    alpha = 0.99; % 平滑因子 (通常需设定一个遗忘因子，取值接近1)
    
    for n = L:N
        % 1. 构建当前时刻的输入向量 (倒序排列)
        xx = x(n:-1:n-L+1);
        
        % 2. 滤波计算
        y(n) = w' * xx;
        e(n) = d(n) - y(n);
        
        % 3. 更新辅助变量 beta(n) 和 p(n)
        beta = alpha * beta + (1 - alpha) * abs(e(n));
        p = alpha * p + (1 - alpha) * (e(n)^2);
        
        % 4. 计算变步长 mu(n)      
        % 计算输入信号功率 x(n)'x(n)

        x_power = xx' * xx;      
        % 计算分子

        numerator = beta * (1 - tanh(beta * abs(e(n))));        
        % 计算分母 (加入 eps 防止除零)

        denominator = x_power + beta * p + eps;        
        mu = numerator / denominator;    
        
        % 5. 权系数更新
        w = w + mu * e(n) * xx;
    end
end