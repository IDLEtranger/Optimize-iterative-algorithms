function [rosenbrock_t, beale_t, goldstein_t] = result_output(iterative_method, epsilon, x0_rosenbrock, x0_beale, x0_goldstein)
    % input: iterative_method——迭代方法
    %        epsilon——允许误差 default: 0.001
    %        x0——迭代初始点 defalut: (0, 0)
    % 为参数赋默认值
    if nargin < 2 || isempty(epsilon)
        epsilon = 0.001;          
    end
    if nargin < 3 || isempty(x0_rosenbrock)
        x0_rosenbrock = [0; 0];
    end
    if nargin < 4 || isempty(x0_beale)
        x0_beale = [4; 0.6];
    end
    if nargin < 5 || isempty(x0_goldstein)
        x0_goldstein = [1; -3];
        % goldstein price函数具有多个局部极小点，需要选取合适初始点才可迭代到全局最优点
    end
    % 迭代计算 
    [rosenbrock_t.x_bar, rosenbrock_t.fmin, rosenbrock_t.data] = iterative_method(@rosenbrock, x0_rosenbrock, epsilon);
    [beale_t.x_bar, beale_t.fmin, beale_t.data] = iterative_method(@beale, x0_beale, epsilon);
    [goldstein_t.x_bar, goldstein_t.fmin, goldstein_t.data] = iterative_method(@goldstein_price, x0_goldstein, epsilon);
    % plot prepare
    x = -5:0.1:5;
    y = -5:0.1:5;
    [X, Y] = meshgrid(x, y);
    
    %%%%%%%%%%%% rosenbrock %%%%%%%%%%%%
    % 绘制函数图像
    figure('Name', strcat('Rosenbrock_', func2str(iterative_method)));
    Z_rosenbrock = zeros(size(X));
    for ir_1 = 1:size(X, 1)
        for jr_1 = 1:size(X, 2)
            Z_rosenbrock(ir_1, jr_1) = rosenbrock([X(ir_1, jr_1); Y(ir_1, jr_1)]);
        end
    end
    surf(X, Y, Z_rosenbrock)
    % 绘制收敛路径
    hold on
    x_r = zeros(1, rosenbrock_t.data.k);
    y_r = zeros(1, rosenbrock_t.data.k);
    z_r = zeros(1, rosenbrock_t.data.k);
    for ir_2 = 1:rosenbrock_t.data.k
        x_r(ir_2) = rosenbrock_t.data.x{ir_2}(1);
        y_r(ir_2) = rosenbrock_t.data.x{ir_2}(2);
        z_r(ir_2) = rosenbrock([x_r(ir_2); y_r(ir_2)]);       
        % 计算散点的颜色
        scatterColor = [1 1-ir_2/rosenbrock_t.data.k 0];
        scatter3(x_r(ir_2), y_r(ir_2), z_r(ir_2), 15, 'filled', 'MarkerFaceColor', scatterColor);        
        if ir_2 > 1
            plot3([x_r(ir_2-1), x_r(ir_2)], [y_r(ir_2-1), y_r(ir_2)], [z_r(ir_2-1), z_r(ir_2)], 'Color', scatterColor, 'LineWidth', 2);
        end
    end
    xlabel("x"); ylabel("y"); zlabel("z");
    % 输出迭代过程表格
    indices = 1:rosenbrock_t.data.k;
    xk = cell2mat(rosenbrock_t.data.x(1:rosenbrock_t.data.k));
    xk_x = xk(1, :);  % x坐标
    xk_y = xk(2, :);  % y坐标
    rosenbrock_t.iteration_process = array2table([indices' z_r' xk_x' xk_y'], 'VariableNames', {'Index', 'Value', 'x', 'y'});

    %%%%%%%%%%%% beale %%%%%%%%%%%%                                         
    figure('Name', strcat('Beale_', func2str(iterative_method)));
    Z_beale = zeros(size(X));
    for ir_1 = 1:size(X, 1)
        for jr_1 = 1:size(X, 2)
            Z_beale(ir_1, jr_1) = beale([X(ir_1, jr_1); Y(ir_1, jr_1)]);
        end
    end
    surf(X, Y, Z_beale)
    % 绘制收敛路径
    hold on
    x_r = zeros(1, beale_t.data.k);
    y_r = zeros(1, beale_t.data.k);
    z_r = zeros(1, beale_t.data.k);
    for ir_2 = 1:beale_t.data.k
        x_r(ir_2) = beale_t.data.x{ir_2}(1);
        y_r(ir_2) = beale_t.data.x{ir_2}(2);
        z_r(ir_2) = beale([x_r(ir_2); y_r(ir_2)]);       
        % 计算散点的颜色
        scatterColor = [1 1-ir_2/beale_t.data.k 0];
        scatter3(x_r(ir_2), y_r(ir_2), z_r(ir_2), 15, 'filled', 'MarkerFaceColor', scatterColor);        
        if ir_2 > 1
            plot3([x_r(ir_2-1), x_r(ir_2)], [y_r(ir_2-1), y_r(ir_2)], [z_r(ir_2-1), z_r(ir_2)], 'Color', scatterColor, 'LineWidth', 2);
        end
    end
    xlabel("x"); ylabel("y"); zlabel("z");
    % 输出迭代过程表格
    indices = 1:beale_t.data.k;
    xk = cell2mat(beale_t.data.x(1:beale_t.data.k));
    xk_x = xk(1, :);  % x坐标
    xk_y = xk(2, :);  % y坐标
    beale_t.iteration_process = array2table([indices' z_r' xk_x' xk_y'], 'VariableNames', {'Index', 'Value', 'x', 'y'});

    %%%%%%%%%%%% goldstein_price %%%%%%%%%%%%
    figure('Name', strcat('Goldstein_price_', func2str(iterative_method)));
    Z_goldstein = zeros(size(X));
    for ir_1 = 1:size(X, 1)
        for jr_1 = 1:size(X, 2)
            Z_goldstein(ir_1, jr_1) = goldstein_price([X(ir_1, jr_1); Y(ir_1, jr_1)]);
        end
    end
    surf(X, Y, Z_goldstein)

    % 绘制收敛路径
    hold on
    x_r = zeros(1, goldstein_t.data.k);
    y_r = zeros(1, goldstein_t.data.k);
    z_r = zeros(1, goldstein_t.data.k);
    for ir_2 = 1:goldstein_t.data.k
        x_r(ir_2) = goldstein_t.data.x{ir_2}(1);
        y_r(ir_2) = goldstein_t.data.x{ir_2}(2);
        z_r(ir_2) = goldstein_price([x_r(ir_2); y_r(ir_2)]);       
        % 计算散点的颜色
        scatterColor = [1 1-ir_2/goldstein_t.data.k 0];
        scatter3(x_r(ir_2), y_r(ir_2), z_r(ir_2), 15, 'filled', 'MarkerFaceColor', scatterColor);        
        if ir_2 > 1
            plot3([x_r(ir_2-1), x_r(ir_2)], [y_r(ir_2-1), y_r(ir_2)], [z_r(ir_2-1), z_r(ir_2)], 'Color', scatterColor, 'LineWidth', 2);
        end
    end
    xlabel("x"); ylabel("y"); zlabel("z");

    % 输出迭代过程表格
    indices = 1:goldstein_t.data.k;
    xk = cell2mat(goldstein_t.data.x(1:goldstein_t.data.k));
    xk_x = xk(1, :);  % x坐标
    xk_y = xk(2, :);  % y坐标
    goldstein_t.iteration_process = array2table([indices' z_r' xk_x' xk_y'], 'VariableNames', {'Index', 'Value', 'x', 'y'});
