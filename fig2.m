%% fig2.m：论文图 2，噪声与滞回、时间演化和节点图样
% 布局: (a)转变曲线，(b)阈值统计，(c)滞后宽度统计，(d)时间轨迹。
clear; clc; close all;
script_dir=fileparts(mfilename('fullpath'));
addpath(script_dir,fullfile(script_dir,'simulations'));

%% ============================================================
%% 初始化参数
%% ============================================================
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;

% Figure 2a 参数
ETA_LIST = [0, 0.9, 1.0];
DYNA.sigma_min  = 0;
DYNA.sigma_max  = 24;
DYNA.sigma_npts = 121;
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

% Figure 2d 参数
DYNA_CD.sigma = 16.0;
scan_mode = 'bwd';
noise_list = [0.00, 0.80];

%% ============================================================
%% 创建主图窗和布局 (2x2 网格)
%% ============================================================
hf = figure('Visible', 'on', 'Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.94, 0.9]);
tl_main = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

plots_dir = fullfile(script_dir, 'fig');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

%% ============================================================
%% (a) 上左 (1,1): Figure 2a - A vs Sigma 曲线
%% ============================================================
ax_a = nexttile(tl_main, 1);
hold on;
colors = lines(length(ETA_LIST));

% 加载并绘制正向扫描 (Forward)
for i = 1:length(ETA_LIST)
    cur_eta = ETA_LIST(i);
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, cur_eta);
    data_file = fullfile(script_dir, 'results', mat_name);
    if exist(data_file, 'file')
        res = load(data_file);
        plot(res.sigma_range, res.A_fwd , '-o', 'Color', colors(i,:), 'LineWidth', 3, ...
            'MarkerSize', 7.5, 'MarkerFaceColor', colors(i,:), ...
            'MarkerIndices', 1:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Fwd, $\\eta = %.1f$', cur_eta));
    end
end

% 加载并绘制反向扫描 (Backward)
for i = 1:length(ETA_LIST)
    cur_eta = ETA_LIST(i);
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, cur_eta);
    data_file = fullfile(script_dir, 'results', mat_name);
    if exist(data_file, 'file')
        res = load(data_file);
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 3, ...
            'MarkerSize', 7.5, ...
            'MarkerIndices', 3:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Bwd, $\\eta = %.1f$', cur_eta));
    end
end

xlabel('$\sigma$', 'FontSize', 27, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$A(\sigma)$', 'FontSize', 27, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 24]);
% ylim([0, 148]);
legend('Location', 'NorthWest', 'FontSize', 27, 'NumColumns', 1, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
grid on;
set(ax_a, 'FontSize', 27, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
% 添加面板标签 (a)
annotation(hf, 'textbox', [0.005, 0.933, 0.05, 0.05], 'String', '(a)', ...
    'FontSize', 27, 'FontWeight', 'bold', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

%% ============================================================
%% (b) 上右 (1,2): Figure 1b - 三相图
%% ============================================================
ax_b = nexttile(tl_main, 2);
hold(ax_b, 'on');
DYNA_ETA = DYNA;
DYNA_ETA.noise = 0.01;
DYNA_ETA.T_END = 500;
DYNA_ETA.steps = 2;
DYNA_ETA.init_perturb = 0.1;
DYNA_ETA.sigma_min = 0;
DYNA_ETA.sigma_max = 35;

P_VAL = TOPO_PARAM;
ETA_LIST = linspace(0, 2, 201);
SEED_LIST = 1:5;
FORCE_RERUN = false;
if numel(SEED_LIST) < 2 || numel(unique(SEED_LIST)) ~= numel(SEED_LIST) || ...
        any(SEED_LIST < 0 | SEED_LIST ~= floor(SEED_LIST))
    error('SEED_LIST 至少需要两个互不重复的非负整数种子。');
end

script_dir = fileparts(mfilename('fullpath'));
res_dir = fullfile(script_dir, 'results');
num_eta = numel(ETA_LIST);
[SF_REPEATS, SB_REPEATS] = load_noise_repeats( ...
    script_dir, DYNA_ETA, TOPO_TYPE, P_VAL, ETA_LIST, SEED_LIST, FORCE_RERUN);

eta_range = ETA_LIST;
sf_vec = mean(SF_REPEATS, 1);
sb_vec = mean(SB_REPEATS, 1);
sf_var = var(SF_REPEATS, 0, 1);
sb_var = var(SB_REPEATS, 0, 1);
sf_std = sqrt(sf_var);
sb_std = sqrt(sb_var);

stats_path = fullfile(res_dir, sprintf( ...
    'scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f_threshold_stats.mat', ...
    P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));
save(stats_path, 'P_VAL', 'ETA_LIST', 'SEED_LIST', 'SF_REPEATS', 'SB_REPEATS', ...
    'sf_vec', 'sb_vec', 'sf_var', 'sb_var', 'sf_std', 'sb_std');

max_val = max([sf_vec + sf_std, sb_vec + sb_std]) * 1.15;
min_val = min(sb_vec) * 0.85;
fill_max_vec = ones(size(eta_range)) * max_val;
fill_min_vec = ones(size(eta_range)) * min_val;

% 三个区域由两条均值曲线分界，与原图的定义相同。
fill(ax_b, [eta_range, fliplr(eta_range)], [sb_vec, fill_min_vec], ...
    [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
fill(ax_b, [eta_range, fliplr(eta_range)], [sf_vec, fliplr(sb_vec)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off');
fill(ax_b, [eta_range, fliplr(eta_range)], [fill_max_vec, fliplr(sf_vec)], ...
    [1.0 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');

% 均值 ± 1 个样本标准差；标准差与纵轴 sigma 的单位相同。
forward_color = [0 0.447 0.741];
backward_color = [0.85 0.325 0.098];
fill(ax_b, [eta_range, fliplr(eta_range)], ...
    [sf_vec + sf_std, fliplr(sf_vec - sf_std)], forward_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.18, 'HandleVisibility', 'off');
fill(ax_b, [eta_range, fliplr(eta_range)], ...
    [sb_vec + sb_std, fliplr(sb_vec - sb_std)], backward_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.18, 'HandleVisibility', 'off');

plot(ax_b, eta_range, sf_vec, '--', 'Color', forward_color, 'LineWidth', 3.75, ...
    'DisplayName', 'Forward $\sigma_c^f$');
plot(ax_b, eta_range, sb_vec, '--', 'Color', backward_color, 'LineWidth', 3.75, ...
    'DisplayName', 'Backward $\sigma_c^b$');

t1 = text(ax_b, 0.8, 10, 'Uniform', ...
    'HorizontalAlignment', 'center', 'FontSize', 27, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
t2 = text(ax_b, 0.23, 17, 'Hysteresis', ...
    'HorizontalAlignment', 'center', 'FontSize', 27, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
t3 = text(ax_b, 1.6, 14, 'Turing', ...
    'HorizontalAlignment', 'center', 'FontSize', 27, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

xlabel(ax_b, '$\eta$', 'FontSize', 21, 'Interpreter', 'latex');
ylabel(ax_b, '$\sigma$', 'FontSize', 27, 'Interpreter', 'latex');
grid(ax_b, 'off'); box(ax_b, 'on');
xlim(ax_b, [0, 2]); ylim(ax_b, [min_val, max_val]);
set(ax_b, 'FontSize', 27, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

legend(ax_b, 'Location', 'northeast', 'FontSize', 27, 'Interpreter', 'latex', 'FontName', 'Times New Roman');

grid(ax_b, 'off');
box(ax_b, 'on');
set(ax_b, 'FontSize', 27, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
% 添加面板标签 (b)
annotation(hf, 'textbox', [0.485, 0.917, 0.05, 0.05], 'String', '(b)', ...
    'FontSize', 27, 'FontWeight', 'bold', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

%% ============================================================
%% (c) 下左 (2,1): 四种连接概率的平均滞后宽度及样本标准差
%% ============================================================
ax_c = nexttile(tl_main, 3);
hold(ax_c, 'on');
P_LIST = [0.01, 0.03, 0.05, 0.07];
width_colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556];
width_mean = zeros(numel(P_LIST), num_eta);
width_std = zeros(numel(P_LIST), num_eta);
for p_idx = 1:numel(P_LIST)
    p_value = P_LIST(p_idx);
    if abs(p_value-P_VAL) < 1e-10
        sf_samples = SF_REPEATS;
        sb_samples = SB_REPEATS;
    else
        [sf_samples, sb_samples] = load_noise_repeats( ...
            script_dir, DYNA_ETA, TOPO_TYPE, p_value, ETA_LIST, SEED_LIST, FORCE_RERUN);
    end
    % 先逐次配对计算宽度，再计算均值和样本标准差，保留阈值间的协方差。
    delta_samples = sf_samples - sb_samples;
    width_mean(p_idx,:) = mean(delta_samples, 1);
    width_std(p_idx,:) = std(delta_samples, 0, 1);
    color = width_colors(p_idx,:);
    fill(ax_c, [ETA_LIST, fliplr(ETA_LIST)], ...
        [width_mean(p_idx,:)+width_std(p_idx,:), ...
        fliplr(width_mean(p_idx,:)-width_std(p_idx,:))], color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax_c, ETA_LIST, width_mean(p_idx,:), 'Color', color, 'LineWidth', 3, ...
        'DisplayName', sprintf('$p=%.2f$', p_value));
end
xlabel(ax_c, '$\eta$', 'FontSize', 27, 'Interpreter', 'latex');
ylabel(ax_c, '$\langle\Delta\sigma\rangle$', 'FontSize', 27, 'Interpreter', 'latex');
xlim(ax_c, [0 2]);
box(ax_c, 'on');
set(ax_c, 'FontSize', 27, 'FontName', 'Times New Roman', 'Layer', 'top', ...
    'TickLabelInterpreter', 'latex');
legend(ax_c, 'Location', 'northeast', 'FontSize', 23, 'Interpreter', 'latex', ...
    'NumColumns', 2, 'FontName', 'Times New Roman');
annotation(hf, 'textbox', [0.005, 0.443, 0.05, 0.05], 'String', '(c)', ...
    'FontSize', 27, 'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', 'Times New Roman');
save(fullfile(res_dir, 'fig2_hysteresis_width_stats.mat'), ...
    'P_LIST', 'ETA_LIST', 'SEED_LIST', 'width_mean', 'width_std');

%% ============================================================
%% (d) 下右 (2,2): eta=0.0 和 0.8 的时间轨迹
%% ============================================================
U_data = cell(2, 1);
t_data = cell(2, 1);

% 加载轨迹数据
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
for i = 1:2
    cur_noise = noise_list(i);
    mat_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s%.1f_n%.2f_%s_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA_CD.sigma, cur_noise, scan_mode);
    file_path = fullfile(results_dir, mat_name);

    if exist(file_path, 'file')
        data = load(file_path, 't', 'Y');
        t_data{i} = data.t;
        U_data{i} = data.Y(:, 1:DYNA.N);
    else
        warning('未找到数据文件: %s', mat_name);
    end
end

% 自动探测特定节点 (突变轨线)
idx_switch = [];
if ~isempty(U_data{2})
    u2 = U_data{2};
    u2_end = u2(end, :);
    u2_early = u2(round(size(u2, 1)*0.2), :);
    idx_switch = find(u2_early > 4 & u2_end < 2, 1);
    if isempty(idx_switch)
        [~, idx_switch] = max(max(u2) - min(u2));
    end
end

% 绘制轨迹
c_u = [0.0 0.447 0.741];
% 嵌套布局
tl_c = tiledlayout(tl_main, 2, 1, 'Padding', 'none', 'TileSpacing', 'compact');
tl_c.Layout.Tile = 4;
for i = 1:2
    ax_c_sub = nexttile(tl_c);
    hold(ax_c_sub, 'on');

    if ~isempty(U_data{i})
        t_cur = t_data{i};
        u_cur = U_data{i};
        u_end = u_cur(end, :);
        [~, idx_up] = max(u_end);
        [~, idx_down] = min(u_end);

        plot(ax_c_sub, t_cur, u_cur(:, 1:5:end), 'Color', [0.85 0.85 0.85], 'LineWidth', 0.75);
        plot(ax_c_sub, t_cur, u_cur(:, idx_up), 'Color', c_u, 'LineWidth', 2.25);
        plot(ax_c_sub, t_cur, u_cur(:, idx_down), 'Color', c_u*0.6, 'LineWidth', 2.25);
        if i == 2 && ~isempty(idx_switch)
            plot(ax_c_sub, t_cur, u_cur(:, idx_switch), 'Color', [0.9 0.4 0], 'LineWidth', 3.0);
        end
    end

    ylabel(ax_c_sub, '$u(t)$', 'FontSize', 27, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    xlim(ax_c_sub, [0, 500]);
    grid(ax_c_sub, 'on');
    set(ax_c_sub, 'FontSize', 27, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');

    if ~isempty(U_data{i})
        text(ax_c_sub, 0.95, 0.8, sprintf('$\\eta = %.1f$', noise_list(i)), ...
            'Units', 'normalized', 'FontSize', 27, 'Interpreter', 'latex', ...
            'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 2, 'HorizontalAlignment', 'right', 'FontName', 'Times New Roman');
    end

    if i == 2
        xlabel(ax_c_sub, '$t$', 'FontSize', 27, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    end
end

% 添加时间轨迹面板标签 (d)
annotation(hf, 'textbox', [0.485, 0.440, 0.05, 0.05], 'String', '(d)', ...
    'FontSize', 27, 'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

% 在布局完成后添加近零宽度放大图。
drawnow;
pos_c = ax_c.Position;
ax_inset = axes(hf, 'Position', [pos_c(1)+0.54*pos_c(3), ...
    pos_c(2)+0.36*pos_c(4), 0.40*pos_c(3), 0.31*pos_c(4)]);
hold(ax_inset, 'on');
for p_idx = 1:numel(P_LIST)
    color = width_colors(p_idx,:);
    fill(ax_inset, [ETA_LIST, fliplr(ETA_LIST)], ...
        [width_mean(p_idx,:)+width_std(p_idx,:), ...
        fliplr(width_mean(p_idx,:)-width_std(p_idx,:))], color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(ax_inset, ETA_LIST, width_mean(p_idx,:), 'Color', color, 'LineWidth', 1.5);
end
xlim(ax_inset, [1.5 1.6]);
ylim(ax_inset, [-0.12 0.14]);
xticks(ax_inset, [1.5 1.55 1.6]);
yticks(ax_inset, [-0.1 0 0.1]);
set(ax_inset, 'FontSize', 16, 'FontName', 'Times New Roman', 'Box', 'on', ...
    'Layer', 'top', 'TickLabelInterpreter', 'latex');

out_fig = fullfile(plots_dir, 'fig2.eps');
exportgraphics(hf, out_fig, 'ContentType', 'vector');
fprintf('[SUCCESS] 综合图表已导出: %s\n', out_fig);


function [sf_samples, sb_samples] = load_noise_repeats(script_dir, dyna, topology, p, eta_list, seeds, force_rerun)
% 每种概率按指定种子复用缓存；只有缺失文件才补算，统计次数由 SEED_LIST 决定。
res_dir = fullfile(script_dir, 'results');
sf_samples = zeros(numel(seeds), numel(eta_list));
sb_samples = zeros(size(sf_samples));
for seed_idx = 1:numel(seeds)
    seed = seeds(seed_idx);
    data_path = fullfile(res_dir, sprintf( ...
        'scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f_seed%.0f.mat', ...
        p, lower(topology), dyna.N, dyna.K, dyna.alpha, dyna.beta, seed));
    if force_rerun || ~exist(data_path, 'file')
        fprintf('[SIM] 补算 p=%.3f、seed=%d。\n', p, seed);
        dyna_seed = dyna;
        dyna_seed.noise_seed = seed;
        sigma_eta_connectivity_er(dyna_seed, topology, p, eta_list);
    else
        fprintf('[CACHE] 复用 p=%.3f、seed=%d。\n', p, seed);
    end
    data = load(data_path, 'p_val', 'ETA_LIST', 'sf_vec', 'sb_vec');
    if ~isfield(data, 'p_val') || abs(data.p_val-p) > 1e-10 || ...
            ~isfield(data, 'ETA_LIST') || ~isequal(data.ETA_LIST(:), eta_list(:)) || ...
            ~isfield(data, 'sf_vec') || numel(data.sf_vec) ~= numel(eta_list) || ...
            ~isfield(data, 'sb_vec') || numel(data.sb_vec) ~= numel(eta_list) || ...
            any(~isfinite(data.sf_vec(:))) || any(~isfinite(data.sb_vec(:)))
        error('p=%.3f、seed=%d 的缓存与当前参数不符：%s', p, seed, data_path);
    end
    sf_samples(seed_idx,:) = data.sf_vec(:)';
    sb_samples(seed_idx,:) = data.sb_vec(:)';
end
end
