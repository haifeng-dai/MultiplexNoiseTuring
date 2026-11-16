%% fig3.m
% Combine spectral reference lines and simulation points in a 2x3 layout.
% The deterministic spectral reference predicts forward instability only; it excludes noise.
% Row 1: (beta/alpha, sigma) phase diagrams for selected alpha values.
% Row 2: (alpha, sigma) phase diagrams for selected beta/alpha ratios.
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

%% 1. Config
TOPO_TYPE = 'ER';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
MARKER_STEP_ABC = 5;  % Display every fifth simulation point in the top row.
MARKER_STEP_DEF = 5;  % Display every fifth simulation point in the bottom row.
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
plots_dir = fullfile(script_dir, 'fig');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% --- Row 1 Config (alpha values) ---
SELECTED_ALPHAS = [0.03, 0.05, 0.10];
ALPHA_TOL = 1e-10; % Floating-point tolerance only.
data_path_abc = fullfile(res_dir, sprintf('scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));

% alpha_config for (a, b, c)
alpha_config(1).alpha = 0.03;
alpha_config(1).region_x_norms = [0.7, 0.35, 0.75]; alpha_config(1).region_y_norms = [0.10, 0.3, 0.38];

alpha_config(2).alpha = 0.05;
alpha_config(2).region_x_norms = [0.7, 0.35, 0.6]; alpha_config(2).region_y_norms = [0.10, 0.2, 0.70];

alpha_config(3).alpha = 0.10;
alpha_config(3).region_x_norms = [0.55, 0.7, 0.55]; alpha_config(3).region_y_norms = [0.1, 0.5, 0.92];

% --- Row 2 Config (ratio values) ---
SELECTED_RATIOS = [0.1, 6.0, 10.0];
RATIO_TOL = 1e-10; % Floating-point tolerance only.
data_path_def = fullfile(res_dir, sprintf('scan_alpha_ratio_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));

% ratio_config for (d, e, f)
ratio_config(1).ratio = 0.1;
ratio_config(1).region_x_norms = [0.7, 0.4, 0.6]; ratio_config(1).region_y_norms = [0.08, 0.2, 0.72];

ratio_config(2).ratio = 6.0;
ratio_config(2).region_x_norms = [0.5, 0.83, 0.35]; ratio_config(2).region_y_norms = [0.2, 0.47, 0.75];

ratio_config(3).ratio = 10.0;
ratio_config(3).region_x_norms = [0.8, 0.65, 0.8]; ratio_config(3).region_y_norms = [0.08, 0.5, 0.9];

%% 2. Load Data
if ~exist(data_path_abc, 'file') || ~exist(data_path_def, 'file')
    error('Missing data files. Ensure sigma_beta_alpha.m and sigma_alpha_beta.m have been run.');
end
data_abc = load(data_path_abc);
data_def = load(data_path_def);

%% 3. Plotting
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.55]);

% Helpers

% --- Combined 2x3 Layout ---
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
axs = gobjects(6, 1);

% Row 1
for i = 1:3
    axs(i) = nexttile(t); hold on;
    target_alpha = SELECTED_ALPHAS(i);
    idx_a = find(abs(data_abc.ALPHA_LIST - target_alpha) <= ALPHA_TOL);
    if isempty(idx_a)
        error('fig2:MissingAlpha', 'Requested alpha=%.4g is absent from %s. Available values: %s', ...
            target_alpha, data_path_abc, mat2str(data_abc.ALPHA_LIST));
    elseif numel(idx_a) ~= 1
        error('fig2:AmbiguousAlpha', 'Multiple entries match alpha=%.4g in %s.', target_alpha, data_path_abc);
    end
    actual_alpha = data_abc.ALPHA_LIST(idx_a);
    x_data = data_abc.RATIO_LIST(:)';
    sf_vec = data_abc.SF_Matrix(idx_a, :);
    sb_vec = data_abc.SB_Matrix(idx_a, :);
    sf_sim = sf_vec; sb_sim = sb_vec;
    [sigma_spec, lambda_star] = forward_spectral_threshold( ...
        TOPO_TYPE, N, K, P_VAL_FIXED, actual_alpha, x_data);
    mask = sf_vec < sb_vec;
    if any(mask), [sf_vec(mask), sb_vec(mask)] = deal(sb_vec(mask), sf_vec(mask)); end
    max_val = max([sf_sim, sb_sim, sigma_spec(isfinite(sigma_spec))]) * 1.12; min_val = max(0, min(sb_vec) * 0.88);
    if i == 1, min_val = 8; end
    fill_max_vec = ones(size(x_data)) * max_val; fill_min_vec = ones(size(x_data)) * min_val;
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    % Background regions use simulation thresholds; points retain scan identity.
    visible_idx = find(x_data >= 0 & x_data <= 10);
    marker_idx = unique([visible_idx(1:MARKER_STEP_ABC:end), visible_idx(end)]);
    h_spec = plot(x_data, sigma_spec, '-', 'Color', [0.00 0.447 0.741], ...
        'LineWidth', 1.8, 'DisplayName', 'Theory $\sigma_s^f$');
    h_forward = plot(x_data, sf_sim, 'o', 'LineStyle', 'none', ...
        'Color', [0.00 0.447 0.741], 'MarkerFaceColor', [0.00 0.447 0.741], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 7.5, 'MarkerIndices', marker_idx, ...
        'DisplayName', 'Fwd $\sigma_c^f$');
    h_backward = plot(x_data, sb_sim, 's', 'LineStyle', 'none', ...
        'Color', [0.85 0.325 0.098], 'MarkerFaceColor', [0.85 0.325 0.098], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 7.5, 'MarkerIndices', marker_idx, ...
        'DisplayName', 'Bwd $\sigma_c^b$');
    xlabel(sprintf('$\\beta/\\alpha$ ($\\alpha = %.2f$)', actual_alpha), 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    if i == 1, ylabel('$\sigma$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman'); end
    xlim([0, 10]); ylim([min_val, max_val]);
    set(axs(i), 'FontSize', 21, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
    if i == 1, legend([h_spec, h_forward, h_backward], 'Location', 'NorthEast', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman'); end
    grid off; box on;
    text(0.05, 0.95, sprintf('\\textbf{(%s)}', char('a'+i-1)), 'Units', 'normalized', 'FontSize', 21, 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontName', 'Times New Roman');

    % Map coordinates for annotations
    drawnow;
    xl = xlim(axs(i)); yl = ylim(axs(i));
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));

    config = alpha_config(i);
    % Region text
    text(norm2data_x(config.region_x_norms(1)), min_val + config.region_y_norms(1)*(max_val-min_val), 'Uniform', 'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(2)), min_val + config.region_y_norms(2)*(max_val-min_val), 'Hysteresis', 'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    if i == 1
        turing_x = 3;
        turing_y = 26;
    else
        turing_x = norm2data_x(config.region_x_norms(3));
        turing_y = min_val + config.region_y_norms(3)*(max_val-min_val);
    end
    text(turing_x, turing_y, 'Turing', 'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
end

% Row 2
for i = 1:3
    axs(i+3) = nexttile(t); hold on;
    target_ratio = SELECTED_RATIOS(i);
    idx_r = find(abs(data_def.RATIO_LIST - target_ratio) <= RATIO_TOL);
    if isempty(idx_r)
        error('fig2:MissingRatio', 'Requested beta/alpha=%.4g is absent from %s. Available values: %s', ...
            target_ratio, data_path_def, mat2str(data_def.RATIO_LIST));
    elseif numel(idx_r) ~= 1
        error('fig2:AmbiguousRatio', 'Multiple entries match beta/alpha=%.4g in %s.', target_ratio, data_path_def);
    end
    actual_ratio = data_def.RATIO_LIST(idx_r);
    x_data = data_def.ALPHA_LIST(:)';
    sf_vec = data_def.SF_Matrix(:, idx_r)'; sb_vec = data_def.SB_Matrix(:, idx_r)';
    sf_sim = sf_vec; sb_sim = sb_vec;
    [sigma_spec, lambda_star] = forward_spectral_threshold( ...
        TOPO_TYPE, N, K, P_VAL_FIXED, x_data, actual_ratio);
    max_val = max([sf_sim, sb_sim, sigma_spec(isfinite(sigma_spec))]) * 1.15; min_val = 0;
    if i == 3, min_val = 10; end
    fill_max_vec = ones(size(x_data)) * max_val; fill_min_vec = ones(size(x_data)) * min_val;
    % Fills
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    % Background regions use simulation thresholds; points retain scan identity.
    visible_idx = find(x_data >= 0.01 & x_data <= 0.1);
    marker_idx = unique([visible_idx(1:MARKER_STEP_DEF:end), visible_idx(end)]);
    h_spec = plot(x_data, sigma_spec, '-', 'Color', [0.00 0.447 0.741], ...
        'LineWidth', 1.8, 'DisplayName', 'Theory $\sigma_c^f$');
    h_forward = plot(x_data, sf_sim, 'o', 'LineStyle', 'none', ...
        'Color', [0.00 0.447 0.741], 'MarkerFaceColor', [0.00 0.447 0.741], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 7.5, 'MarkerIndices', marker_idx, ...
        'DisplayName', 'Fwd $\sigma_c^f$');
    h_backward = plot(x_data, sb_sim, 's', 'LineStyle', 'none', ...
        'Color', [0.85 0.325 0.098], 'MarkerFaceColor', [0.85 0.325 0.098], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 7.5, 'MarkerIndices', marker_idx, ...
        'DisplayName', 'Bwd $\sigma_c^b$');
    xlabel(sprintf('$\\alpha$ ($\\beta/\\alpha = %.1f$)', actual_ratio), 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    if i == 1, ylabel('$\sigma$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman'); end
    xlim([0.01, 0.1]); ylim([min_val, max_val]);
    set(axs(i+3), 'FontSize', 21, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
    grid off; box on;
    text(0.05, 0.95, sprintf('\\textbf{(%s)}', char('d'+i-1)), 'Units', 'normalized', 'FontSize', 21, 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontName', 'Times New Roman');

    % Map coordinates for annotations
    drawnow;
    xl = xlim(axs(i+3)); yl = ylim(axs(i+3));
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));

    config = ratio_config(i);
    % Region text
    text(norm2data_x(config.region_x_norms(1)), min_val + config.region_y_norms(1)*(max_val-min_val), 'Uniform', 'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(2)), min_val + config.region_y_norms(2)*(max_val-min_val), 'Hysteresis', 'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(3)), min_val + config.region_y_norms(3)*(max_val-min_val), 'Turing', 'HorizontalAlignment', 'center', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
end

% Export
out_path = fullfile(plots_dir, 'fig3.eps');
exportgraphics(h, out_path, 'ContentType', 'vector');
