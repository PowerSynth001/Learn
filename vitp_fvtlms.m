function [y, e, L_history] = vitp_fvtlms(d, x, max_L, L_init)
% 输入参数:
% d: 期望信号
% x: 输入信号
% max_L: 滤波器最大允许阶数
% L_init: 初始阶数
%
% 输出参数:
% y: 滤波器输出
% e: 误差信号
% L_history: 阶数变化历史

    N = length(x);
    w = zeros(max_L, 1);       % 权系数初始化
    L_frac = L_init;           % 分数阶数
    L_curr = floor(L_frac);    % 当前整数阶数
    
    y = zeros(N, 1);
    e = zeros(N, 1);
    L_history = zeros(N, 1);
    
    % --- 算法参数初始化 
    Delta = 2;                 % 片段长度差 (用于计算片段误差)
    gamma = 0.002;             % 泄漏因子 (用于阶数收敛)
    beta = 0.99;               % 平滑因子 
    lambda = 0.05;             % 常数xi
    C = 2;                     % 限定函数常数
    alpha_para = 0.01;         % 迭代参数 alpha (初始值)
    mu_w = 0.01;               % 权值更新步长 (固定或变步长均可)

    for n = max_L:N
        % 1. 获取当前整数阶数
        L_curr = floor(L_frac);
        % 限制阶数范围
        L_curr = max(Delta + 1, min(max_L, L_curr));
        
        % 2. 计算完整误差 (Full Error)
        x_full = x(n:-1:n-L_curr+1);
        y_full = w(1:L_curr)' * x_full;
        e_full = d(n) - y_full;
        
        % 保存输出
        y(n) = y_full;
        e(n) = e_full;
        L_history(n) = L_frac;
        
        % 3. 计算片段误差 (Segment Error)
       
        L_seg = L_curr - Delta;
        x_seg = x(n:-1:n-L_seg+1);
        y_seg = w(1:L_seg)' * x_seg;
        e_seg = d(n) - y_seg;
        
        % 4. 计算误差平方差 BE(n)
       
        BE = abs(e_full^2 - e_seg^2);
        
        % 5. 更新迭代参数 alpha(n)
        
        delta_alpha = lambda * BE;
               
        % 注意：此处利用平滑处理避免 alpha 剧烈跳变
        alpha_para = beta * alpha_para + (1 - beta) * delta_alpha;
        
        % 6. 更新分数阶数 L_frac
        % 计算阶数梯度项 (这是FVTLMS的通用梯度近似，此处根据文档描述加入限定函数)
        % 梯度通常包含泄漏项和误差相关项
        % 假设 gradient_term 为标准 FVTLMS 的梯度方向
        
        
        % 简化的梯度估计 (对应文中隐含的阶数调整方向)
        % 这里的 gradient_term 近似为: gamma * L(n) - e(n)^2 (或类似能量项)
        % 为了体现"限定函数"，我们将梯度项输入 tanh
        
        % 构造梯度输入 x_grad (示例逻辑，具体梯度形式视FVTLMS基础模型定)
        % 通常阶数减少的方向与泄漏因子 gamma 有关，增加方向与误差有关
        term_gradient = gamma * L_curr - e_full^2; 
        
        % 应用限定函数 Phi(x)
      
        phi_val = C * tanh(term_gradient / C);
        
        % 更新阶数
       
        L_frac = L_frac - alpha_para * phi_val;
        
        % 边界检查
        L_frac = max(Delta + 2, min(max_L, L_frac));
        
        % 7. 更新权系数 (LMS)
        w(1:L_curr) = w(1:L_curr) + mu_w * e_full * x_full;
        % 超出当前阶数的部分置零或保持衰减
        if L_curr < max_L
            w(L_curr+1:end) = 0;
        end
    end
end