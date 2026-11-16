%% fig4.m：计算、绘制并合并论文图 4
% 原两图的计算和绘图分别封装为本文件局部函数，不依赖原绘图脚本。
% 六个坐标轴和色条直接创建于最终图窗，统一两行高度与边界。
% 默认只读缓存；允许计算时仅补缺失缓存/噪声行，不强制重算。
clear; clc; close all;
RUN_SPECTRAL_COMPUTATION=false;
RUN_INTERACTION_COMPUTATION=false;
script_dir=fileparts(mfilename('fullpath'));
addpath(script_dir,fullfile(script_dir,'simulations'));
fig_dir=fullfile(script_dir,'fig');
if ~isfolder(fig_dir),mkdir(fig_dir);end
FIGURE_VISIBLE='on';
fig=figure('Visible',FIGURE_VISIBLE,'Color','w','Units','normalized', ...
    'Position',[0.03 0.05 0.90 0.60]);
set(fig,'Units','pixels','Position',[40 60 1500 640]);
% 一次性定义最终布局；左侧为谱图，右侧为四个参数平面。
row_y=[0.59 0.13];row_height=0.33;
positions=[0.065 row_y(1) 0.42 row_height; ...
    0.065 row_y(2) 0.36 row_height; ...
    0.535 row_y(1) 0.165 row_height; ...
    0.78 row_y(1) 0.165 row_height; ...
    0.535 row_y(2) 0.165 row_height; ...
    0.78 row_y(2) 0.165 row_height];
panel_axes=gobjects(6,1);
for idx=1:6
    panel_axes(idx)=axes('Parent',fig,'Units','normalized', ...
        'Position',positions(idx,:),'PositionConstraint','innerposition');
end
renderSpectrum(panel_axes(1:2),script_dir,RUN_SPECTRAL_COMPUTATION);
renderInteraction(panel_axes(3:6),script_dir,RUN_INTERACTION_COMPUTATION);
% 色条也使用最终位置；恢复坐标轴位置，避免自动色条布局缩轴。
bar_panels=[2 3 4 5 6];bar_x=[0.435 0.71 0.955 0.71 0.955];
for idx=1:numel(bar_panels)
    panel=bar_panels(idx);cb=colorbar(panel_axes(panel));
    cb.Units='normalized';
    cb.Position=[bar_x(idx) positions(panel,2) 0.009 row_height];
    if panel==2
        cb.Label.String='Eigenvector Component';
        set(cb.Label,'FontSize',14,'FontName','Times New Roman','Interpreter','latex');
    end
end
for idx=1:6,panel_axes(idx).Position=positions(idx,:);end
% 合并后显式统一文字字号，避免标签倍率和区域缩放造成差异。
axes_list=findall(fig,'Type','axes');
for idx=1:numel(axes_list)
    ax=axes_list(idx);
    if any(ax==panel_axes(3:6))
        axes_font_size=14;label_font_size=16;
    else
        axes_font_size=21;label_font_size=21;
    end
    set(ax,'FontSize',axes_font_size,'FontName','Times New Roman','TickLabelInterpreter','latex');
    set([ax.XLabel ax.YLabel],'FontSize',label_font_size,'FontName','Times New Roman', ...
        'Interpreter','latex','FontWeight','normal');
end
colorbars=findall(fig,'Type','colorbar');
for idx=1:numel(colorbars)
    cb=colorbars(idx);
    colorbar_font_size=14;
    set(cb,'FontSize',colorbar_font_size,'FontName','Times New Roman','TickLabelInterpreter','latex');
    set(cb.Label,'FontSize',21,'FontName','Times New Roman','Interpreter','latex');
    if strcmp(cb.Label.String,'Eigenvector Component')
        cb.Label.FontSize=14;
        cb.Label.FontSizeMode='manual';
    end
end
drawnow;
% 标号放在 y label 顶部外侧，顶部低于坐标框顶端。
for panel_idx=1:numel(panel_axes)
    letter=char('a'+panel_idx-1);ax=panel_axes(panel_idx);
    ax.YLabel.Units='normalized';label_position=ax.YLabel.Position;
    ax.Units='pixels';axis_position=ax.Position;
    label_x=max(label_position(1)-0.02,(16-axis_position(1))/axis_position(3));
    text(ax,label_x,0.98,['\textbf{(' letter ')}'], ...
        'Units','normalized','FontSize',21,'FontName','Times New Roman', ...
        'Interpreter','latex','HorizontalAlignment','center', ...
        'VerticalAlignment','top','Clipping','off');
end
% 隐藏交互工具栏，避免图窗截图出现省略号按钮。
all_axes=findall(fig,'Type','axes');
for idx=1:numel(all_axes),all_axes(idx).Toolbar.Visible='off';end
output_file=fullfile(fig_dir,'fig4.eps');
% 普通 figure 中只有坐标轴和色条，统一导出整张图。
set(fig,'PaperPositionMode','auto');
print(fig,output_file,'-depsc','-painters');
print(fig,fullfile(fig_dir,'fig4.png'),'-dpng','-r300');
fprintf('[完成] 图 4：%s\n',output_file);
if strcmp(FIGURE_VISIBLE,'off'),close(fig);end

function renderSpectrum(panel_axes,script_dir,run_computation)
%% 1. Load existing spectral results
res_dir = fullfile(script_dir, 'results');

input_file = fullfile(res_dir, ...
    'compare_forward_threshold_spectral_er_N200_K5_p0.030_a0.100.mat');
if ~exist(input_file, 'file')
    assert(run_computation,'缺少谱缓存；请启用 RUN_SPECTRAL_COMPUTATION：%s',input_file);
    computeSpectralReference(script_dir);
end
S = load(input_file);

RATIO_LIST = getRequiredField(S, {'RATIO_LIST', 'ratio_list'});
lambda_star = getRequiredField(S, {'lambda_star', 'LAMBDA_STAR'});
lambda_all = getRequiredField(S, {'lambda_all', 'LAMBDA_ALL'});
phi_star = getRequiredField(S, {'phi_star', 'PHI_STAR'});
RATIO_LIST = RATIO_LIST(:)';
lambda_star = lambda_star(:)';
if isfield(S, 'N'), N = double(S.N); else, N = 200; end
if isfield(S, 'K'), K = double(S.K); else, K = 5; end
if isfield(S, 'ALPHA_FIXED')
    alpha = double(S.ALPHA_FIXED);
elseif isfield(S, 'alpha')
    alpha = double(S.alpha);
else
    alpha = 0.1;
end
if abs(alpha-0.1) > 1e-12
    error('Input alpha=%.12g does not match the required alpha=0.1.', alpha);
end
if numel(lambda_star) ~= numel(RATIO_LIST)
    error('lambda_star and RATIO_LIST must have the same length.');
end

%% 2. Compute the optimal spectral position of the continuous dispersion relation
lambda_upper = (10/3) / alpha;
sigma_of_lambda = @(lambda) (110/3 + 4*alpha.*lambda) ./ ...
    ((alpha.*lambda) .* (10/3 - alpha.*lambda));
lambda_grid = linspace(lambda_upper*1e-8, lambda_upper*(1-1e-8), 100001);
sigma_grid = sigma_of_lambda(lambda_grid);
[~, grid_min_idx] = min(sigma_grid);
bracket_lo = lambda_grid(max(1, grid_min_idx-1));
bracket_hi = lambda_grid(min(numel(lambda_grid), grid_min_idx+1));
lambda_opt = fminbnd(sigma_of_lambda, bracket_lo, bracket_hi);

%% 3. Load the nondegenerate dominant eigenvector at r=4.5
target_ratio = 4.5;
[ratio_error, target_idx] = min(abs(RATIO_LIST-target_ratio));
if ratio_error > 1e-7
    error('r=4.5 was not found in RATIO_LIST.');
end
target_ratio = RATIO_LIST(target_idx);
target_lambda = lambda_star(target_idx);

spectrum = real(getSpectrumAt(lambda_all, target_idx));
spectrum = spectrum(:);
degeneracy_tolerance = max(1e-8, 1e-6*max(1, abs(target_lambda)));
multiplicity = sum(abs(spectrum-target_lambda) < degeneracy_tolerance);
if multiplicity ~= 1
    warning('The dominant eigenvalue at r=%.1f has multiplicity %d; a nondegenerate eigenvalue was expected.', ...
        target_ratio, multiplicity);
end

phi = getEigenvectorAt(phi_star, target_idx);
if numel(phi) ~= N*K || norm(phi) == 0
    error('Invalid phi_star vector at r=%.1f.', target_ratio);
end
phi = real(phi(:));
phi = phi / norm(phi, 2);
[~, imax] = max(abs(phi));
if phi(imax) < 0
    phi = -phi;
end
M = reshape(phi, [N, K])';
IPR_full = sum(abs(phi).^4) / (sum(abs(phi).^2)^2);

%% 4. Plot the combined figure
selected_ratios = [4.3 4.4 4.5 4.6];
selected_indices = zeros(size(selected_ratios));
for j = 1:numel(selected_ratios)
    [ratio_error, selected_indices(j)] = min(abs(RATIO_LIST-selected_ratios(j)));
    if ratio_error > 1e-7
        warning('r=%.1f is absent from RATIO_LIST; its marker is omitted.', selected_ratios(j));
        selected_indices(j) = NaN;
    end
end



% Top panel: reuse the lambda_star curve from dominant_mode_switch_alpha010.
ax1 = panel_axes(1);
ax1.Tag='panel_a';
plot(ax1, RATIO_LIST, lambda_star, '-o', 'LineWidth', 1.5, ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.0 0.447 0.741]);
hold(ax1, 'on');
yline(ax1, lambda_opt, '--', sprintf('$\\lambda_{\\mathrm{opt}}=%.4g$', lambda_opt), ...
    'LineWidth', 1.4, 'LabelHorizontalAlignment', 'left', ...
    'FontSize', 21, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
for j = 1:numel(selected_ratios)
    if isfinite(selected_indices(j))
        idx = selected_indices(j);
        plot(ax1, RATIO_LIST(idx), lambda_star(idx), 'kp', ...
            'MarkerFaceColor', [1.0 0.75 0.0], 'MarkerSize', 10, ...
            'HandleVisibility', 'off');
    end
end
xlim(ax1, [min(RATIO_LIST) max(RATIO_LIST)]);
xlabel(ax1, '$\beta/\alpha$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel(ax1, '$\lambda_*$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
grid(ax1, 'on'); box(ax1, 'on');
set(ax1, 'FontSize', 21, 'FontName', 'Times New Roman', 'TickLabelInterpreter', 'latex');

% Bottom panel: reuse the signed eigenvector heatmap at r=4.5 from the bottom-right panel of com2_3.png.
% Keep the heatmap, colorbar, and colorbar text within the top-panel width.
ax2 = panel_axes(2);
ax2.Tag='panel_b';
imagesc(ax2, 1:N, 1:K, M);
set(ax2, 'YDir', 'normal', 'YTick', 1:K);
signed_map = blueWhiteRedMap(257);
signed_color_limit = 0.4;
clim(ax2, [-signed_color_limit signed_color_limit]);
colormap(ax2, signed_map);
xlabel(ax2, 'Node index $i$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel(ax2, 'Layer $k$', 'FontSize', 21, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
set(ax2, 'FontSize', 21, 'FontName', 'Times New Roman', 'TickLabelInterpreter', 'latex');
% Reserve additional columns for the bottom-panel colorbar and its text.
cb = colorbar(ax2, 'Location', 'eastoutside');
cb.Label.String = 'Eigenvector Component';
set(cb, 'FontSize', 21, 'FontName', 'Times New Roman', 'TickLabelInterpreter', 'latex');
cb.Label.FontSize = 14;
cb.Label.FontName = 'Times New Roman';
cb.Label.Interpreter = 'latex';
box(ax2, 'on');



end

function renderInteraction(panel_axes,script_dir,run_computation)
project_dir = script_dir;
addpath(fullfile(project_dir, 'simulations'));

%% 1. 固定参数与扫描变量
TOPO_TYPE = 'ER';
TOPO_PARAM = 0.03;
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 30;

etaList = 0:0.1:2;
ratioList = 0:0.1:10.0;

%% 2. Reuse cached results and simulate missing eta rows only
numEta = numel(etaList);
numRatio = numel(ratioList);
result_dir = fullfile(project_dir, 'results');
file_stem = sprintf('eta_beta_interaction_%s_N%d_K%d_p%.3f_a%.3f', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha);
result_path = fullfile(result_dir, [file_stem, '.mat']);
assembled_rows = false;
simulation_ready = false;

if exist(result_path, 'file')
    data = load(result_path, 'etaList', 'ratioList', 'sigmaF', 'sigmaB');
    if ~isfield(data, 'etaList') || ~isfield(data, 'ratioList') || ...
            ~isequal(data.etaList(:), etaList(:)) || ...
            ~isequal(data.ratioList(:), ratioList(:)) || ...
            ~isfield(data, 'sigmaF') || ~isfield(data, 'sigmaB') || ...
            ~isequal(size(data.sigmaF), [numEta, numRatio]) || ...
            ~isequal(size(data.sigmaB), [numEta, numRatio]) || ...
            any(~isfinite(data.sigmaF(:))) || any(~isfinite(data.sigmaB(:)))
        error('Invalid or incompatible cached results: %s', result_path);
    end
    sigmaF = data.sigmaF;
    sigmaB = data.sigmaB;
    fprintf('[LOAD] %s\n', result_path);
else
    assembled_rows = true;
    if ~exist(result_dir, 'dir'), mkdir(result_dir); end
    sigmaF = nan(numEta, numRatio);
    sigmaB = nan(numEta, numRatio);
    for etaIdx = 1:numEta
        eta_result_path = fullfile(result_dir, sprintf('%s_eta%.2f.mat', file_stem, etaList(etaIdx)));
        if ~exist(eta_result_path, 'file')
            assert(run_computation,'缺少联合参数缓存；请启用 RUN_INTERACTION_COMPUTATION：%s',eta_result_path);
            % Initialize simulation resources only when a row is missing.
            if ~simulation_ready
                net_path = fullfile(project_dir, 'results', 'topology', upper(TOPO_TYPE), ...
                    sprintf('N%d_p%.3f.mat', DYNA.N, TOPO_PARAM));
                if ~exist(net_path, 'file')
                    error('Missing network required to simulate uncached rows: %s', net_path);
                end
                net_data = load(net_path, 'nets');
                L_intra_lib = net_data.nets;
                y_universal_seed = load_standard_pattern_seed(DYNA);
                adj_inter = ones(DYNA.K) - eye(DYNA.K);
                DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;
                if isempty(gcp('nocreate'))
                    parpool('Threads', 64);
                end
                simulation_ready = true;
            end
            fprintf('[SIM] Computing missing eta=%.3f (%d/%d).\n', etaList(etaIdx), etaIdx, numEta);
            sigmaFRow = zeros(1, numRatio);
            sigmaBRow = zeros(1, numRatio);
            parfor ratioIdx = 1:numRatio
                dyna_local = DYNA;
                dyna_local.noise = etaList(etaIdx);
                dyna_local.beta = ratioList(ratioIdx) * dyna_local.alpha;
                dyna_local.L_intra = L_intra_lib;
                dyna_local.y_seed = y_universal_seed;
                thresholds = find_thresholds(TOPO_TYPE, TOPO_PARAM, dyna_local);
                sigmaFRow(ratioIdx) = thresholds(1);
                sigmaBRow(ratioIdx) = thresholds(2);
            end
            if any(~isfinite(sigmaFRow)) || any(~isfinite(sigmaBRow))
                error('Simulation returned invalid thresholds for eta=%.3f.', etaList(etaIdx));
            end
            save(eta_result_path, 'etaList', 'ratioList', 'etaIdx', 'sigmaFRow', 'sigmaBRow', ...
                'TOPO_TYPE', 'TOPO_PARAM', 'DYNA');
            fprintf('[SAVE] %s\n', eta_result_path);
        end
        data = load(eta_result_path, 'etaList', 'ratioList', 'etaIdx', 'sigmaFRow', 'sigmaBRow');
        if ~isfield(data, 'etaList') || ~isfield(data, 'ratioList') || ...
                ~isfield(data, 'etaIdx') || ~isscalar(data.etaIdx) || ...
                data.etaIdx < 1 || data.etaIdx > numel(data.etaList) || ...
                data.etaIdx ~= floor(data.etaIdx) || ...
                abs(data.etaList(data.etaIdx)-etaList(etaIdx)) > 1e-10 || ...
                ~isequal(data.ratioList(:), ratioList(:)) || ...
                ~isfield(data, 'sigmaFRow') || ~isfield(data, 'sigmaBRow') || ...
                ~isequal(size(data.sigmaFRow), [1, numRatio]) || ...
                ~isequal(size(data.sigmaBRow), [1, numRatio]) || ...
                any(~isfinite(data.sigmaFRow(:))) || any(~isfinite(data.sigmaBRow(:)))
            error('Invalid or incompatible cached row: %s', eta_result_path);
        end
        sigmaF(etaIdx, :) = data.sigmaFRow;
        sigmaB(etaIdx, :) = data.sigmaBRow;
        fprintf('[LOAD] eta=%.3f: %s\n', etaList(etaIdx), eta_result_path);
    end
end

deltaSigma = sigmaF - sigmaB;

etaZeroIdx = find(abs(etaList) < eps, 1);
ratioZeroIdx = find(abs(ratioList) < eps, 1);
if isempty(etaZeroIdx) || isempty(ratioZeroIdx)
    error('interactionIndex 计算要求 etaList 和 ratioList 均包含 0。');
end
interactionIndex = sigmaF ...
    - repmat(sigmaF(etaZeroIdx, :), numEta, 1) ...
    - repmat(sigmaF(:, ratioZeroIdx), 1, numRatio) ...
    + sigmaF(etaZeroIdx, ratioZeroIdx);

assert(isequal(size(sigmaF), [numEta, numRatio]));
assert(isequal(size(sigmaB), [numEta, numRatio]));
assert(isequal(size(deltaSigma), [numEta, numRatio]));
assert(isequal(size(interactionIndex), [numEta, numRatio]));
% Save a complete cache after all rows have been loaded or computed.
if assembled_rows
    save(result_path, 'etaList', 'ratioList', 'sigmaF', 'sigmaB', ...
        'deltaSigma', 'interactionIndex', 'TOPO_TYPE', 'TOPO_PARAM', 'DYNA');
    fprintf('[SAVE] Complete results: %s\n', result_path);
end

%% 5. 2 × 2 参数平面图（仅显示实际扫描点，不进行插值或平滑）
[ratioGrid, etaGrid] = meshgrid(ratioList, etaList);
ax1 = panel_axes(1);
ax1.Tag='panel_c';
scatter(ax1, ratioGrid(:), etaGrid(:), 150, sigmaF(:), 's', 'filled', 'MarkerEdgeColor', 'none');
xlabel(ax1, '$\beta/\alpha$', 'Interpreter', 'latex'); ylabel(ax1, '$\eta$', 'Interpreter', 'latex');
colorbar(ax1); grid(ax1, 'off'); box(ax1, 'on');

ax2 = panel_axes(2);
ax2.Tag='panel_d';
scatter(ax2, ratioGrid(:), etaGrid(:), 150, sigmaB(:), 's', 'filled', 'MarkerEdgeColor', 'none');
xlabel(ax2, '$\beta/\alpha$', 'Interpreter', 'latex'); ylabel(ax2, '$\eta$', 'Interpreter', 'latex');
colorbar(ax2); grid(ax2, 'off'); box(ax2, 'on');

ax3 = panel_axes(3);
ax3.Tag='panel_e';
scatter(ax3, ratioGrid(:), etaGrid(:), 150, deltaSigma(:), 's', 'filled', 'MarkerEdgeColor', 'none');
xlabel(ax3, '$\beta/\alpha$', 'Interpreter', 'latex'); ylabel(ax3, '$\eta$', 'Interpreter', 'latex');
colorbar(ax3); grid(ax3, 'off'); box(ax3, 'on');

ax4 = panel_axes(4);
ax4.Tag='panel_f';
scatter(ax4, ratioGrid(:), etaGrid(:), 150, interactionIndex(:), 's', 'filled', ...
    'MarkerEdgeColor', 'none');
xlabel(ax4, '$\beta/\alpha$', 'Interpreter', 'latex'); ylabel(ax4, '$\eta$', 'Interpreter', 'latex');
colorbar(ax4); grid(ax4, 'off'); box(ax4, 'on');
interactionLimit = max(abs(interactionIndex(:)));
if interactionLimit == 0, interactionLimit = 1; end
clim(ax4, [-interactionLimit, interactionLimit]);
halfMap = 128;
blueWhiteRed = [linspace(0.15, 1, halfMap)', linspace(0.30, 1, halfMap)', ones(halfMap, 1); ...
    ones(halfMap, 1), linspace(1, 0.20, halfMap)', linspace(1, 0.15, halfMap)'];
colormap(ax4, blueWhiteRed);

xPad = 0.02 * max(ratioList(end) - ratioList(1), 1);
yPad = 0.02 * max(etaList(end) - etaList(1), 1);
xLimits = [ratioList(1) - xPad, ratioList(end) + xPad];
yLimits = [etaList(1) - yPad, etaList(end) + yPad];
set([ax1, ax2, ax3, ax4], 'FontSize', 21, 'FontName', 'Times New Roman', 'TickLabelInterpreter', 'latex', 'Layer', 'top', ...
    'XLim', xLimits, 'YLim', yLimits);
set([ax1, ax2, ax3, ax4], 'XTick', 0:2:10, 'YTick', 0:0.5:2, 'XTickLabelRotation', 0);
% Match com2_3_1.m typography for labels and colorbars.
axes_list = [ax1, ax2, ax3, ax4];
for axIdx = 1:numel(axes_list)
    ax = axes_list(axIdx);
    set([ax.XLabel, ax.YLabel], 'FontSize', 21, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    cb = colorbar(ax);
    set(cb, 'FontSize', 21, 'FontName', 'Times New Roman', 'TickLabelInterpreter', 'latex');
end
yLabels = [ax1.YLabel, ax2.YLabel, ax3.YLabel, ax4.YLabel];
set(yLabels, 'Units', 'normalized');
for labelIdx = 1:numel(yLabels)
    labelPos = yLabels(labelIdx).Position;
    labelPos(1) = -0.14;
    yLabels(labelIdx).Position = labelPos;
end


end

function computeSpectralReference(script_dir)
%% 1. Fig. 2 参数与输出配置
addpath(script_dir, fullfile(script_dir, 'simulations'));

TOPO_TYPE = 'ER';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
ALPHA_FIXED = 0.1;
ETA_REFERENCE = 0.01;  % Fig. 2 仿真配置中的噪声参数；谱计算不使用它

ALPHA_TOL = 1e-10;
NONZERO_TOL = 1e-10;
ROOT_TOL = 1e-6;
COARSE_STEP = 0.1;
SIGMA_INITIAL_MAX = 100;
SIGMA_MAX_LIMIT = 2000;

res_dir = fullfile(script_dir, 'results');
fig_dir = fullfile(script_dir, 'fig');

% Fig. 2 的 beta/alpha 扫描结果由 sigma_beta_alpha.m 写入。
threshold_file = fullfile(res_dir, sprintf( ...
    'scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));
topology_file = fullfile(res_dir, 'topology', upper(TOPO_TYPE), ...
    sprintf('N%d_p%.3f.mat', N, P_VAL_FIXED));

if ~exist(threshold_file, 'file')
    error('找不到 Fig. 2 forward 阈值文件：%s', threshold_file);
end
if ~exist(topology_file, 'file')
    error('找不到 Fig. 2 网络文件：%s', topology_file);
end

%% 2. 读取 Fig. 2 forward 阈值和 beta/alpha 横坐标
threshold_data = load(threshold_file);
required_threshold_fields = {'SF_Matrix', 'ALPHA_LIST'};
for field_idx = 1:numel(required_threshold_fields)
    field_name = required_threshold_fields{field_idx};
    if ~isfield(threshold_data, field_name)
        error('Fig. 2 阈值文件缺少字段：%s', field_name);
    end
end

if isfield(threshold_data, 'RATIO_LIST')
    RATIO_LIST = threshold_data.RATIO_LIST(:)';
elseif isfield(threshold_data, 'BETA_LIST')
    if ~isfield(threshold_data, 'ALPHA_LIST')
        error('阈值文件只有 BETA_LIST，但缺少 ALPHA_LIST，无法恢复 beta/alpha。');
    end
    RATIO_LIST = threshold_data.BETA_LIST(:)' / ALPHA_FIXED;
else
    error('Fig. 2 阈值文件缺少 RATIO_LIST 或 BETA_LIST。');
end

alpha_list = threshold_data.ALPHA_LIST(:)';
[alpha_error, alpha_idx] = min(abs(alpha_list - ALPHA_FIXED));
if alpha_error > ALPHA_TOL
    error('Fig. 2 结果中没有 alpha=%.8g，最近值为 %.8g。', ...
        ALPHA_FIXED, alpha_list(alpha_idx));
end

sigma_c_sim = threshold_data.SF_Matrix(alpha_idx, :);
sigma_c_sim = sigma_c_sim(:)';
if numel(sigma_c_sim) ~= numel(RATIO_LIST)
    error('RATIO_LIST 与 SF_Matrix 的列数不一致。');
end

if isfield(threshold_data, 'P_VAL_FIXED')
    if abs(threshold_data.P_VAL_FIXED - P_VAL_FIXED) > 1e-12
        error('阈值文件的 p=%.8g，与当前配置 p=%.8g 不一致。', ...
            threshold_data.P_VAL_FIXED, P_VAL_FIXED);
    end
end

fprintf('Fig. 2 forward 阈值文件：%s\n', threshold_file);
fprintf('Fig. 2 网络 realization：%s\n', topology_file);
fprintf('参数：N=%d，K=%d，alpha=%.8g，eta=%.8g，p=%.8g\n', ...
    N, K, ALPHA_FIXED, ETA_REFERENCE, P_VAL_FIXED);
fprintf('读取到 %d 个 beta/alpha 点，alpha 行索引为 %d。\n', ...
    numel(RATIO_LIST), alpha_idx);

%% 3. 读取与 Fig. 2 完全相同的 ER 网络 realization
topology_data = load(topology_file, 'nets');
if ~isfield(topology_data, 'nets') || numel(topology_data.nets) < K
    error('网络文件中的层数不足：需要 %d 层。', K);
end

L_intra = cell(K, 1);
for layer_idx = 1:K
    Lk = topology_data.nets{layer_idx};
    % 当前项目的 gen_er/gen_ws/gen_ba/gen_sf_gamma 保存的就是 D-A。
    % 若误传邻接矩阵，行和通常不为零，这里直接报错，避免静默得到错误谱。
    if max(abs(full(sum(Lk, 2)))) > 1e-8 || any(diag(Lk) < -1e-12)
        error('nets{%d} 不是拉普拉斯矩阵。若源文件保存邻接矩阵，应先构造 D-A。', layer_idx);
    end
    if norm(full(Lk - Lk'), 'fro') > 1e-8
        error('nets{%d} 不是对称拉普拉斯矩阵，无法使用当前实谱分析。', layer_idx);
    end
    L_intra{layer_idx} = sparse(Lk);
end

adj_inter = ones(K) - eye(K);
L_inter = diag(sum(adj_inter, 2)) - adj_inter;
L_intra_all = blkdiag(L_intra{:});
L_inter_big = kron(L_inter, speye(N));

%% 4. 计算每个 beta/alpha 的谱阈值
% Mimura-Murray 稳态 (u*,v*)=(5,10) 处的 Jacobian。
J = [10 / 3, -5; 10, -4];
lambda_all = cell(size(RATIO_LIST));
lambda_nonzero_all = cell(size(RATIO_LIST));
sigma_c_spec = NaN(size(RATIO_LIST));
lambda_star = NaN(size(RATIO_LIST));
q_star = NaN(size(RATIO_LIST));
mu_star = NaN(size(RATIO_LIST));
mu_at_sim_threshold = NaN(size(RATIO_LIST));
phi_star = cell(size(RATIO_LIST));

for ratio_idx = 1:numel(RATIO_LIST)
    ratio = RATIO_LIST(ratio_idx);

    % 论文记号：L = L^L + (beta/alpha) L^I。
    L_super = L_intra_all + ratio * L_inter_big;
    [eigenvectors, eigenvalues_matrix] = eig(full(L_super));
    [lambda_sorted, order] = sort(real(diag(eigenvalues_matrix)), 'ascend');
    eigenvectors = eigenvectors(:, order);

    if lambda_sorted(1) < -1e-8 || any(lambda_sorted < -1e-8)
        error('ratio=%.8g 的 supra-Laplacian 出现明显负特征值。', ratio);
    end
    if abs(lambda_sorted(1)) > 1e-8
        warning('ratio=%.8g 的最小特征值不是数值零：lambda_1=%.8g。', ...
            ratio, lambda_sorted(1));
    end

    nonzero_mask = lambda_sorted > NONZERO_TOL;
    lambda_nonzero = lambda_sorted(nonzero_mask);
    if isempty(lambda_nonzero)
        warning('ratio=%.8g 没有非零空间模态，谱阈值设为 NaN。', ratio);
        lambda_all{ratio_idx} = lambda_sorted;
        lambda_nonzero_all{ratio_idx} = lambda_nonzero;
        continue;
    end

    lambda_all{ratio_idx} = lambda_sorted;
    lambda_nonzero_all{ratio_idx} = lambda_nonzero;

    finite_sim_thresholds = sigma_c_sim(isfinite(sigma_c_sim));
    if isempty(finite_sim_thresholds)
        sigma_upper = SIGMA_INITIAL_MAX;
    else
        sigma_upper = max(SIGMA_INITIAL_MAX, max(finite_sim_thresholds) * 1.25);
    end

    [sigma_root, mu_root, root_found] = find_spectral_threshold( ...
        lambda_nonzero, ALPHA_FIXED, J, 0, sigma_upper, ...
        SIGMA_MAX_LIMIT, COARSE_STEP, ROOT_TOL);

    if ~root_found
        warning(['ratio=%.8g 在 sigma <= %.8g 范围内没有从负到正的线性谱失稳，' ...
            'sigma_c_spec 设为 NaN。'], ratio, SIGMA_MAX_LIMIT);
        continue;
    end

    sigma_c_spec(ratio_idx) = sigma_root;
    mu_star(ratio_idx) = mu_root;

    [~, q_local, mu_modes] = spectral_growth( ...
        sigma_root, lambda_nonzero, ALPHA_FIXED, J);
    lambda_star(ratio_idx) = lambda_nonzero(q_local);
    nonzero_indices = find(nonzero_mask);
    q_star(ratio_idx) = nonzero_indices(q_local);
    phi_star{ratio_idx} = eigenvectors(:, q_star(ratio_idx));

    if isfinite(sigma_c_sim(ratio_idx))
        [mu_at_sim_threshold(ratio_idx), ~, ~] = spectral_growth( ...
            sigma_c_sim(ratio_idx), lambda_nonzero, ALPHA_FIXED, J);
    end

    % 记录根附近的符号检查，确认是从稳定侧穿过 mu=0。
    delta_sigma = max(1e-4, 10 * ROOT_TOL);
    [mu_below, ~, ~] = spectral_growth( ...
        max(0, sigma_root - delta_sigma), lambda_nonzero, ALPHA_FIXED, J);
    [mu_above, ~, ~] = spectral_growth( ...
        sigma_root + delta_sigma, lambda_nonzero, ALPHA_FIXED, J);
    if ~(mu_below < 0 && mu_above > 0)
        warning(['ratio=%.8g 的根附近符号检查未满足 mu_below<0<mu_above：' ...
            'below=%.8g，root=%.8g，above=%.8g。'], ...
            ratio, mu_below, mu_root, mu_above);
    end
end

%% 5. 汇总结果并保存 MAT、CSV
difference = sigma_c_sim - sigma_c_spec;
relative_difference = difference ./ sigma_c_sim;

threshold_table = table(RATIO_LIST(:), sigma_c_sim(:), sigma_c_spec(:), ...
    difference(:), relative_difference(:), lambda_star(:), q_star(:), ...
    mu_star(:), mu_at_sim_threshold(:), ...
    'VariableNames', {'ratio', 'sigma_c_sim', 'sigma_c_spec', 'difference', ...
    'relative_difference', 'lambda_star', 'q_star', 'mu_star', ...
    'mu_at_sim_threshold'});

base_name = sprintf('compare_forward_threshold_spectral_%s_N%d_K%d_p%.3f_a%.3f', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED, ALPHA_FIXED);
out_mat = fullfile(res_dir, [base_name, '.mat']);
out_csv = fullfile(res_dir, [base_name, '.csv']);
save(out_mat, 'threshold_table', 'RATIO_LIST', 'sigma_c_sim', 'sigma_c_spec', ...
    'difference', 'relative_difference', 'lambda_all', 'lambda_nonzero_all', ...
    'lambda_star', 'q_star', 'mu_star', 'mu_at_sim_threshold', 'phi_star', ...
    'J', 'ALPHA_FIXED', 'ETA_REFERENCE', 'P_VAL_FIXED', ...
    'NONZERO_TOL', 'ROOT_TOL', 'threshold_file', 'topology_file', '-v7');
writetable(threshold_table, out_csv);


end

function value = getRequiredField(S, candidates)
% Read a required field using candidate names.
for i = 1:numel(candidates)
    if isfield(S, candidates{i})
        value = S.(candidates{i});
        return;
    end
end
error('MAT file is missing a required field. Accepted names: %s', strjoin(candidates, ', '));
end

function spectrum = getSpectrumAt(lambda_all, idx)
% Support cell arrays and spectra stored by ratio in columns or rows.
if iscell(lambda_all)
    spectrum = lambda_all{idx};
elseif size(lambda_all, 2) >= idx
    spectrum = lambda_all(:, idx);
elseif size(lambda_all, 1) >= idx
    spectrum = lambda_all(idx, :)';
else
    error('No spectrum at index %d in lambda_all.', idx);
end
end

function vector = getEigenvectorAt(phi_all, idx)
% Support cell arrays and eigenvectors stored by ratio in columns.
if iscell(phi_all)
    vector = phi_all{idx};
elseif size(phi_all, 2) >= idx
    vector = phi_all(:, idx);
else
    error('No eigenvector at index %d in phi_star.', idx);
end
end

function map = blueWhiteRedMap(n_colors)
% Build a blue-white-red diverging colormap with a light cream color at zero.
half_count = floor((n_colors+1)/2);
red_count = n_colors-half_count+1;
neutral_color = [1.00 0.985 0.95];
blue_to_neutral = [linspace(0, neutral_color(1), half_count)', ...
    linspace(0, neutral_color(2), half_count)', ...
    linspace(1, neutral_color(3), half_count)'];
neutral_to_red = [linspace(neutral_color(1), 1, red_count)', ...
    linspace(neutral_color(2), 0, red_count)', ...
    linspace(neutral_color(3), 0, red_count)'];
map = [blue_to_neutral; neutral_to_red(2:end, :)];
end

function [mu_max, q_star, mu_modes] = spectral_growth(sigma, lambda_nonzero, alpha, J)
% 计算所有非零空间模态的增长率，并返回最危险模态。
a = J(1, 1) - alpha * lambda_nonzero;
b = J(1, 2);
c = J(2, 1);
d = J(2, 2) - alpha * sigma * lambda_nonzero;
trace_modes = a + d;
discriminant = (a - d).^2 + 4 * b * c;
root_discriminant = sqrt(complex(discriminant, 0));
eigenvalue_plus = (trace_modes + root_discriminant) / 2;
eigenvalue_minus = (trace_modes - root_discriminant) / 2;
mu_modes = max(real(eigenvalue_plus), real(eigenvalue_minus));
[mu_max, q_star] = max(mu_modes);
end

function [sigma_root, mu_root, root_found] = find_spectral_threshold( ...
    lambda_nonzero, alpha, J, sigma_min, sigma_initial_max, ...
    sigma_max_limit, coarse_step, root_tol)
% 粗扫描寻找最早的负到正区间，再用二分法精确求根。
sigma_root = NaN;
mu_root = NaN;
root_found = false;
sigma_upper = max(sigma_initial_max, sigma_min + coarse_step);

while true
    sigma_grid = sigma_min:coarse_step:sigma_upper;
    mu_grid = zeros(size(sigma_grid));
    for sigma_idx = 1:numel(sigma_grid)
        [mu_grid(sigma_idx), ~, ~] = spectral_growth( ...
            sigma_grid(sigma_idx), lambda_nonzero, alpha, J);
    end

    crossing_idx = find(mu_grid(1:end-1) <= 0 & mu_grid(2:end) >= 0, 1, 'first');
    if ~isempty(crossing_idx)
        lo = sigma_grid(crossing_idx);
        hi = sigma_grid(crossing_idx + 1);
        mu_lo = mu_grid(crossing_idx);

        for iteration = 1:100
            mid = 0.5 * (lo + hi);
            [mu_mid, ~, ~] = spectral_growth(mid, lambda_nonzero, alpha, J);
            if mu_mid >= 0
                hi = mid;
            else
                lo = mid;
                mu_lo = mu_mid;
            end
            if hi - lo <= root_tol
                break;
            end
        end

        sigma_root = 0.5 * (lo + hi);
        [mu_root, ~, ~] = spectral_growth(sigma_root, lambda_nonzero, alpha, J);
        root_found = isfinite(mu_root) && mu_lo <= 0;
        return;
    end

    if sigma_upper >= sigma_max_limit
        return;
    end
    sigma_upper = min(2 * sigma_upper, sigma_max_limit);
end
end
