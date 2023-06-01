function [x_bar, fmin, data] = steepest_descent(func, x0, epsilon)
    % output：x_bar——最优点 fmin——最优值 
    % input：func——目标函数 x0——初始点 epsilon——允许误差
    
    % 为迭代中的变量预先分配内存
    memallo = 10000;
    d = cell(1, memallo);
    x = cell(1, memallo);
    lambda = zeros(1, memallo);
    % 为x赋初值x0，迭代索引k从1开始递增
    x{1} = x0;
    k = 1;
    % 初始值计算
    [~, grad_tmp] = func(x{k});
    while norm(grad_tmp) > epsilon && k <= memallo
        % 不满足迭代终止条件进入循环, 迭代次数超出预先分配单元大小迭代也停止
        % 更新搜索方向
        d{k} = -grad_tmp;
        % 沿搜索方向做一维搜索
        func_handle = @(lambda) func(x{k} + lambda*d{k});
        lambda(k) = fminsearch(func_handle, 0);
        % 更新序列
        x{k+1} = x{k} + lambda(k)*d{k};
        % 一次迭代完成，迭代索引k自增1
        k = k + 1;
        % 更新梯度
        [~, grad_tmp] = func(x{k});
    end
    if k > memallo
        error("WARNING: 最速下降法迭代次数过多");
    end
    % 将预先分配内存的多余部分删去（删去末尾的空元素）
    % cellfun(@(x) ~isempty(x), Mycell) 用于判断元素是否非空，随后与"1"比较
    lastNonZeroIndex = find(cellfun(@(e) ~isempty(e), d), true, 'last'); 
    d = d(1:lastNonZeroIndex);
    lastNonZeroIndex = find(cellfun(@(e) ~isempty(e), x), true, 'last'); 
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