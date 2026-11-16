%% fig6.m: 以 1x3 布局导出五层纯 WS、纯 BA 与 2WS+2ER+1BA 网络的 (beta, sigma) 三相图
% 曲线与字体对齐既有图；三个网络案例在同一图窗纵向排列。
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, 'simulations'));

RUN_SIMULATION = true;   % 缺少新配置数据时自动补种子并扫描
FORCE_RERUN    = false;  % 是否强制重新计算已存在的数据

%% 1. 全局物理与动力学配置 (完全对齐图 2c)
N = 200;                % 节点数
K = 5;                  % 五层网络；层间全连通
DYNA.N = N;
DYNA.K = K;
DYNA.alpha = 0.10;      % 对齐图 2c 的 alpha = 0.10
DYNA.noise = 0.1;       % 与图5统一的噪声强度
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 80;

% 横轴参数：层间耦合比 ratio = beta / alpha
RATIO_LIST = 0:0.25:10.0; % 41 个采样点
num_r = length(RATIO_LIST);
WS_INITIAL_NEIGHBORS = 6;
WS_REWIRE = 0.01;
BA_M = 3;
ER_P = 0.030;
NETWORK_CONFIG = struct('ws_initial_neighbors',WS_INITIAL_NEIGHBORS, ...
    'ws_rewire',WS_REWIRE,'ba_m',BA_M,'er_p',ER_P);
cache_tag = sprintf('noise%.3f_wsK%d_pr%.2f_bam%d_erP%.3f', ...
    DYNA.noise,WS_INITIAL_NEIGHBORS,WS_REWIRE,BA_M,ER_P);

res_dir = fullfile(script_dir, 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
test_plots_dir = fullfile(script_dir, 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end

%% 2. 定义三种对比网络案例与各自独立的绘图标注坐标配置
% ========== 为每个网络独立配置标注文字和箭头坐标 (0~1 归一化) ==========
% 案例 1: Pure WS (5层全 WS，初始邻居数6，pr=0.01)
CASES(1).name = 'WS';
CASES(1).label = 'Pure WS Multiplex (5 layers)';
CASES(1).tag = 'pure_ws';
CASES(1).suffix = 'ws_K5';
CONFIG(1).region_x_norms = [0.55, 0.70, 0.55]; % Uniform, Hysteresis, Pattern
CONFIG(1).region_y_norms = [0.10, 0.50, 0.92];

% 案例 2: Pure BA (5层全 BA，m=3)
CASES(2).name = 'BA';
CASES(2).label = 'Pure BA Multiplex (5 layers)';
CASES(2).tag = 'pure_ba';
CASES(2).suffix = 'ba_K5';
CONFIG(2).region_x_norms = [0.55, 0.70, 0.55];
CONFIG(2).region_y_norms = [0.10, 0.50, 0.92];

% 案例 3: 第1、2层为 WS，第3、4层为 ER，第5层为 BA
CASES(3).name = 'HETERO';
CASES(3).label = '2WS+2ER+1BA Heterogeneous Multiplex';
CASES(3).tag = 'hetero_2ws_2er_1ba';
CASES(3).suffix = 'hetero_2ws_2er_1ba_K5';
CONFIG(3).region_x_norms = [0.55, 0.70, 0.55];
CONFIG(3).region_y_norms = [0.10, 0.50, 0.92];

%% 3. 并行仿真区 (独立计算与落盘)
scan_data_paths = cell(1, 3);
scan_missing = false(1, 3);
for c_idx = 1:3
    mat_name = sprintf('scan_beta_phase_%s_N%d_K%d_a%.3f_%s.mat', ...
        CASES(c_idx).tag,N,K,DYNA.alpha,cache_tag);
    scan_data_paths{c_idx} = fullfile(res_dir, mat_name);
    scan_missing(c_idx) = ~exist(scan_data_paths{c_idx}, 'file');
end

% 预加载三种底层网络
tmp_ba = load(fullfile(res_dir,'topology','BA',sprintf('N%d_m%d.mat',N,BA_M)),'nets');
tmp_ws = load(fullfile(res_dir,'topology','WS', ...
    sprintf('N%d_K%d_pr%.2f.mat',N,WS_INITIAL_NEIGHBORS,WS_REWIRE)),'nets');
tmp_er = load(fullfile(res_dir,'topology','ER',sprintf('N%d_p%.3f.mat',N,ER_P)),'nets');
if ~all([numel(tmp_ba.nets), numel(tmp_ws.nets), numel(tmp_er.nets)] >= K)
    error('fig6:InsufficientLayers', ...
        'BA、WS、ER 拓扑文件都必须至少包含 %d 层。', K);
end

NET_SETS(1).L = tmp_ws.nets(1:K); % Pure WS
NET_SETS(2).L = tmp_ba.nets(1:K); % Pure BA
NET_SETS(3).L = {tmp_ws.nets{1}, tmp_ws.nets{2}, ...
    tmp_er.nets{1}, tmp_er.nets{2}, tmp_ba.nets{1}}; % 2WS+2ER+1BA

% 准备层间拉普拉斯算子
adj_inter = ones(K) - eye(K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;


% 绘图模式不需要斑图种子。仅当确实要计算缺失/强制重算的数据时准备种子与网络。
need_simulation = RUN_SIMULATION && (FORCE_RERUN || any(scan_missing));
if need_simulation
    seed_file = fullfile(res_dir, sprintf( ...
        'evolution_er_N%d_K%d_p%.3f_a0.050_b0.005_s100.0_n0.00_fwd_results.mat',N,K,ER_P));
    if ~exist(seed_file, 'file')
        fprintf('[SEED] 缺少 N=%d、K=%d 的标准种子；按标准 ER 配置调用 pattern_evolution 补齐。\n', N, K);
        seed_cfg = struct('N', N, 'K', K, 'alpha', 0.05, 'beta', 0.005, ...
            'sigma', 100, 'noise', 0, 'T_END', 500, 'dt', 0.001, ...
            'steps', 2, 'init_perturb', 0.1);
        previous_rng = rng;
        restore_rng = onCleanup(@() rng(previous_rng));
        rng(0, 'twister');
        pattern_evolution(seed_cfg,'ER',ER_P,false);
        clear restore_rng;
        if ~exist(seed_file, 'file')
            error('fig6:SeedGenerationFailed', ...
                'pattern_evolution 未生成预期的标准种子：%s', seed_file);
        end
    end

    % 读取并校验严格匹配 N、K 的种子；不以其他尺寸的状态替代。
    y_universal_seed = load_standard_pattern_seed(DYNA);

    for c_idx = 1:3
        case_info = CASES(c_idx);
        data_path = scan_data_paths{c_idx};

        if FORCE_RERUN || scan_missing(c_idx)
            fprintf('\n============================================================\n');
            fprintf('[SIM] 正在计算: %s (任务 %d/3)...\n', case_info.label, c_idx);
            fprintf('============================================================\n');

            if isempty(gcp('nocreate')), parpool('Threads', 64); end

            L_current = NET_SETS(c_idx).L;
            sf_vec = zeros(1, num_r);
            sb_vec = zeros(1, num_r);
            simStart = tic;

            parfor r_i = 1:num_r
                cur_ratio = RATIO_LIST(r_i);
                dyna_local = DYNA;
                dyna_local.beta = dyna_local.alpha * cur_ratio;
                dyna_local.L_intra = L_current;
                dyna_local.y_seed = y_universal_seed;

                res = find_thresholds('ER', 0.030, dyna_local);
                sf_vec(r_i) = res(1);
                sb_vec(r_i) = res(2);
                fprintf('  [%s] ratio=%.2f 完成 -> sf=%.2f, sb=%.2f\n', ...
                    case_info.name, cur_ratio, res(1), res(2));
            end

            fprintf('[DONE] %s 计算完成，耗时: %.2f 秒。\n', ...
                case_info.label, toc(simStart));
            save(data_path,'sf_vec','sb_vec','RATIO_LIST','DYNA','case_info','NETWORK_CONFIG');
        end
    end
end

%% 4. 绘图与导出：三个案例共用 1x3 布局
fprintf('\n[VIS] 开始绘制 1x3 对比相图...\n');
% 一行子图使用两行图窗的一半高度，保持单个子图尺度接近。
h = figure('Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.275]);
layout = tiledlayout(h,1,3,'TileSpacing','compact','Padding','compact');
axs = gobjects(3,1);

for c_idx = 1:3
    case_info = CASES(c_idx);
    cfg = CONFIG(c_idx);
    mat_name = sprintf('scan_beta_phase_%s_N%d_K%d_a%.3f_%s.mat', ...
        case_info.tag,N,K,DYNA.alpha,cache_tag);
    data_path = fullfile(res_dir, mat_name);

    if ~exist(data_path, 'file')
        warning('未找到数据文件: %s，跳过该图绘制。', mat_name);
        continue;
    end

    data = load(data_path);
    x_data = data.RATIO_LIST;
    sf_vec = data.sf_vec;
    sb_vec = data.sb_vec;

    assert(data.DYNA.N==N && data.DYNA.K==K && ...
        data.DYNA.alpha==DYNA.alpha && data.DYNA.noise==DYNA.noise && ...
        isfield(data,'NETWORK_CONFIG') && isequal(data.NETWORK_CONFIG,NETWORK_CONFIG), ...
        '已有数值数据与当前噪声或网络参数不一致：%s',data_path);

    % 避免填充时出现负宽度
    swap_mask = sf_vec < sb_vec;
    if any(swap_mask)
        tmp = sf_vec(swap_mask);
        sf_vec(swap_mask) = sb_vec(swap_mask);
        sb_vec(swap_mask) = tmp;
    end

    axs(c_idx) = nexttile(layout, c_idx);
    ax = axs(c_idx);
    hold(ax, 'on');

    max_val = max(sf_vec) * 1.12;
    min_val = max(0, min(sb_vec) * 0.88);
    fill_max_vec = ones(size(x_data)) * max_val;
    fill_min_vec = ones(size(x_data)) * min_val;

    % 1. 下层: Uniform (浅青)
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], ...
        [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');

    % 2. 中层: Hysteresis (浅灰)
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], ...
        [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');

    % 3. 上层: Pattern (浅桃粉)
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], ...
        [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');

    % 4. 临界边界虚线
    plot(x_data, sf_vec, '--', 'Color', [0.00 0.447 0.741], 'LineWidth', 1.8, ...
        'DisplayName', 'Forward Critical \sigma_c^f');
    plot(x_data, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.8, ...
        'DisplayName', 'Backward Critical \sigma_c^b');

    % 5. 坐标轴修饰
    xlabel('$\beta/\alpha$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    ylabel('$\sigma$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    box on; grid off;
    xlim([0, 10]);
    ylim([min_val, max_val]);
    set(ax, 'FontSize', 21, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');

    % 6. 区域文字坐标映射与子图标签 (完全对齐 fig3.m)
    text(0.05, 0.95, sprintf('\\textbf{(%s)}', char('a'+c_idx-1)), ...
        'Units', 'normalized', 'FontSize', 21, 'Interpreter', 'latex', ...
        'VerticalAlignment', 'top', 'FontName', 'Times New Roman');
    drawnow;
    xl = xlim(ax);
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));

    % 区域文字 - 使用归一化映射
    t1 = text(norm2data_x(cfg.region_x_norms(1)), min_val + cfg.region_y_norms(1)*(max_val-min_val), 'Uniform', ...
        'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    t2 = text(norm2data_x(cfg.region_x_norms(2)), min_val + cfg.region_y_norms(2)*(max_val-min_val), 'Hysteresis', ...
        'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    t3 = text(norm2data_x(cfg.region_x_norms(3)), min_val + cfg.region_y_norms(3)*(max_val-min_val), 'Turing', ...
        'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');

    uistack([t1, t2, t3], 'top');

end

% 三个案例统一导出。
% test_out_img = fullfile(test_plots_dir,'fig6.png');
% exportgraphics(h,test_out_img,'Resolution',300);
% fprintf('[DONE] 1x3 对比相图已保存：%s\n',test_out_img);
plots_dir=fullfile(script_dir,'fig');
if ~isfolder(plots_dir),mkdir(plots_dir);end
manuscript_eps=fullfile(plots_dir,'fig6.eps');
exportgraphics(h,manuscript_eps,'ContentType','vector');
fprintf('[手稿输出] %s\n',manuscript_eps);

fprintf('\n所有相图处理与导出完毕！\n');
