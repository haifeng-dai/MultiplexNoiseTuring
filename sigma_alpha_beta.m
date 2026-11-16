function sigma_alpha_beta(DYNA, TOPO_TYPE, P_VAL_FIXED, ALPHA_LIST, RATIO_LIST, RUN_SIMULATION, OUTPUT_PATH)
%% sigma_alpha_beta.m: 分析层内扩散(alpha)对滞回阈值(sigma)的影响，并对比不同层间耦合强度(beta)
script_dir=fileparts(mfilename('fullpath'));
addpath(script_dir,fullfile(script_dir,'simulations'));

% 数据保存路径
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
% OUTPUT_PATH 可选；新脚本可指定独立保存路径。
if nargin < 7 || isempty(OUTPUT_PATH)
    data_path = fullfile(res_dir, sprintf('scan_alpha_ratio_%s_N%d_K%d_p%.3f.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, P_VAL_FIXED));
else
    data_path = OUTPUT_PATH;
    output_dir = fileparts(data_path);
    if ~isempty(output_dir) && ~exist(output_dir, 'dir'), mkdir(output_dir); end
end

%% 2. 核心并行计算区 [64核优化版]
num_a = length(ALPHA_LIST);
num_r = length(RATIO_LIST);

if RUN_SIMULATION
    % === 并行池初始化：充分利用所有64个核 ===
    if isempty(gcp('nocreate'))
        parpool('Threads', 64);
        fprintf('[POOL] 已初始化 64 线程并行池\n');
    else
        pool = gcp;
        fprintf('[POOL] 使用现有并行池 (%d 个worker)\n', pool.NumWorkers);
    end

    [R_idx, A_idx] = meshgrid(1:num_r, 1:num_a);
    R_idx = R_idx(:);
    A_idx = A_idx(:);
    total_tasks = length(R_idx);

    Sigma_F_Results = zeros(total_tasks, 1);
    Sigma_B_Results = zeros(total_tasks, 1);

    % --- A. 进度监控器初始化 ---
    q = parallel.pool.DataQueue;
    hWait = waitbar(0, '正在初始化...', 'Name', '并行任务进度');
    cleanupObj = onCleanup(@() delete(hWait));
    tStart_sim = tic;
    afterEach(q, @(~) nUpdateProgress(total_tasks, hWait, tStart_sim));

    % --- B. 构造层间拉普拉斯算子 ---
    adj_inter = (ones(DYNA.K) - eye(DYNA.K));
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    % --- C. 加载全局种子 (保证一致性) ---
    y_universal_seed = load_standard_pattern_seed(DYNA);

    % --- D. 加载固定拓扑网络 ---
    net_path = fullfile(script_dir, 'results', 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, P_VAL_FIXED));
    tmp_net = load(net_path);
    L_intra_lib = tmp_net.nets;

    % 1. 粗略扫描 (Rough Scan) - 用于预测引导
    fprintf('\n[SIM] 执行粗略扫描以提取预测线 (Rough Scan)... \n');
    n_rough_a = 8;
    idx_rough_a = round(linspace(1, num_a, n_rough_a));
    alpha_rough = ALPHA_LIST(idx_rough_a);

    [R_rough, A_rough_idxs] = meshgrid(1:num_r, 1:n_rough_a);
    R_rough = R_rough(:);
    A_rough_idxs = A_rough_idxs(:);
    n_rough_tasks = length(R_rough);

    sf_rough_vec = zeros(n_rough_tasks, 1);
    sb_rough_vec = zeros(n_rough_tasks, 1);

    parfor task_id = 1:n_rough_tasks
        c = DYNA;
        cur_alpha = alpha_rough(A_rough_idxs(task_id));
        c.alpha = cur_alpha;
        c.beta  = cur_alpha * RATIO_LIST(R_rough(task_id));
        c.L_intra = L_intra_lib;
        c.y_seed  = y_universal_seed;
        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, c);
        sf_rough_vec(task_id) = res(1);
        sb_rough_vec(task_id) = res(2);
    end

    sf_rough = reshape(sf_rough_vec, [n_rough_a, num_r]);
    sb_rough = reshape(sb_rough_vec, [n_rough_a, num_r]);

    % 2. 预测插值
    sigma_pred_f = zeros(num_a, num_r);
    sigma_pred_b = zeros(num_a, num_r);
    for r_i = 1:num_r
        sigma_pred_f(:, r_i) = interp1(alpha_rough, sf_rough(:, r_i), ALPHA_LIST, 'pchip', DYNA.sigma_max);
        sigma_pred_b(:, r_i) = interp1(alpha_rough, sb_rough(:, r_i), ALPHA_LIST, 'pchip', DYNA.sigma_max);
    end

    fprintf('\n[SIM] 开始局域精细并行扫描 (Guided Fine Scan): %d 组任务... \n', total_tasks);
    simStart = tic;

    parfor t = 1:total_tasks
        dyna_local = DYNA;
        cur_alpha = ALPHA_LIST(A_idx(t));
        dyna_local.alpha = cur_alpha;
        dyna_local.beta  = cur_alpha * RATIO_LIST(R_idx(t));

        dyna_local.L_intra = L_intra_lib;
        dyna_local.y_seed  = y_universal_seed;

        % 利用预测值缩小搜索区间
        pred_center = max(sigma_pred_f(A_idx(t), R_idx(t)), sigma_pred_b(A_idx(t), R_idx(t)));
        if pred_center > DYNA.sigma_max - 5
            dyna_local.sigma_min = 0;
            dyna_local.sigma_max = DYNA.sigma_max;
        else
            dyna_local.sigma_min = max(0, pred_center - 10);
            dyna_local.sigma_max = min(DYNA.sigma_max, pred_center + 10);
        end

        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, dyna_local);
        Sigma_F_Results(t) = res(1);
        Sigma_B_Results(t) = res(2);
        send(q, t);
    end

    simTime = toc(simStart);
    fprintf('[DONE] 仿真完成！总耗时: %.2f 秒。\n', simTime);

    SF_Matrix = reshape(Sigma_F_Results, [num_a, num_r]);
    SB_Matrix = reshape(Sigma_B_Results, [num_a, num_r]);
    save(data_path, 'SF_Matrix', 'SB_Matrix', 'ALPHA_LIST', 'RATIO_LIST', 'P_VAL_FIXED');
    fprintf('[SAVE] 结果已保存到: %s\n', data_path);
end

% ========================================================
    function nUpdateProgress(total, hWait, tStart)
        persistent count
        if isempty(count), count = 0; end
        count = count + 1;
        progress = count / total;
        waitbar(progress, hWait, sprintf('处理进度: %d/%d (%.0f%%) | 耗时: %.1f 秒', ...
            count, total, progress*100, toc(tStart)));
    end

end
