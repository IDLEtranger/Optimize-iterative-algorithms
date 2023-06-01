function [x_bar, fmin, data] = BFGS(func, x0, epsilon)
    % output：x_bar——最优点 fmin——最优值 
    % input：func——目标函数 x0——初始点 epsilon——允许误差
    
    % 为迭代中需要记录的变量预先分配内存
    memallo = 10000;
    d = cell(1, memallo);
    x = cell(1, memallo);
    lambda = zeros(1, memallo);
    % 为x赋初值x0，迭代索引k从1开始递增
    x{1} = x0;
    % 梯度初始值计算
    [~, grad_tmp] = func(x{1});
    % H g为迭代中间变量无需记录，迭代后便丢弃上一次数据
    Hk_BGFS = eye(length(x{1})); % 此处k=1
    gk = grad_tmp; % 此处k=1
    k = 1;
    
    while true
        % 更新搜索方向
        d{k} = -Hk_BGFS*gk;
        % 沿搜索方向做一维搜索
        func_handle = @(lambda) func(x{k} + lambda*d{k});
        lambda(k) = fminsearch(func_handle, 0);
        % 更新序列
        x{k+1} = x{k} + lambda(k)*d{k};
        % 循环条件判断
        if ~(norm(grad_tmp) > epsilon && k <= memallo)
            break;
        end
        if mod(k, length(x{1})) % 如果距离上次重置已经迭代两次没有跳出循环，则以x{k+1}作为初始点重新开始迭代
            Hk_BGFS = eye(length(x{1})); % 新的初始值
            [~, grad_tmp] = func(x{k+1});
            gk = grad_tmp; % 新的初始值
            k = k + 1;
            continue;
        end
        % 更新梯度
        [~, grad_tmp] = func(x{k+1});
        % 计算迭代中间变量
        gk_next = grad_tmp;
        pk = x{k+1} - x{k};
        qk = gk_next - gk;
        Hk_BGFS_next = Hk_BGFS + (1 + (qk'*Hk_BGFS*qk)/(pk'*qk))*(pk*pk')/(pk'*qk) - (pk*qk'*Hk_BGFS + Hk_BGFS*qk*pk')/(pk'*qk);
        % 一次迭代完成，k自增，同时中间变量完成更新
        k = k + 1;
        Hk_BGFS = Hk_BGFS_next;
        gk = gk_next;
    end
    if k > memallo
        error("WARNING: BFGS拟牛顿法迭代次数过多");
    end
    % 将预先分配内存的多余部分删去（删去末尾的空元素）
    % cellfun(@(x) ~isempty(x), Mycell) 用于判断元素是否非空，随后与"1"比较
    lastNonZeroIndex = find(cellfun(@(e) ~isempty(e), d), 1, 'last'); 
    d = d(1:lastNonZeroIndex);
    lastNonZeroIndex = find(cellfun(@(e) ~isempty(e), x), 1, 'last'); 
    x = x(1:lastNonZeroIndex);
    lastNonZeroIndex = find(lambda, 1, 'last');
    lambda = lambda(1:lastNonZeroIndex);
    % 将中间结果保存在结构体变量中
    data.k = k; % 迭代次数
    data.x = x; % 生成序列
    data.d = d; % 历次搜索方向
    data.lambda = lambda; % 历次步长
    % 函数最优值与最优点
    x_bar = x{k};
    fmin = func(x{k});