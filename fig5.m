%% fig5.m：跨拓扑转变阈值与临界谱点比较（2x3）
% (d) 度与节点数；(e) 层数与层间耦合；(f) 各拓扑枢纽与临界模式关联。
% 只读取下载的网络；每个点为一个 multiplex 实现，不是误差条。
clear; clc; close all;

%% 1. 参数配置
RUN_COMPUTATION = false;       % false：只读取完整谱缓存
RUN_THRESHOLD_COMPUTATION = true;  % 仅补缺失/未完成阈值；已有完成任务只读
RUN_HUB_PROFILE_COMPUTATION = true; % 仅补算缺失的枢纽权重缓存
THRESHOLD_WORKERS = 10;         % 当前缺少 3 个 BA m=1 配置
THRESHOLD_TASK_INDICES = [];  % 空值：全部；可设 1:2 验证少量配置
THRESHOLD_NOISE = 0.1;        % a、b 使用逐边独立乘性噪声
THRESHOLD_NOISE_SEED = 1;     % 固定单个噪声实现，不是统计均值
FIGURE_VISIBLE = 'on';       % 'off'：后台导出图片
ALPHA = 0.05;
COUPLING_RATIO = 0.1;        % beta/alpha
N_VALUES = [100 500];
DEGREE_VALUES = 2:20;
CONNECTIVITY_K = 5;           % a、d 图固定层数、全连通层间耦合
LAYER_VALUES = 3:10;
LAYER_N = 200;               % b、e 图固定节点数与平均层内度
LAYER_DEGREE = 6;
INTERLAYER_TYPES = {'degree2_ring','complete_graph'};
WS_REWIRE = 0.01;            % 标准库 WS 重连概率
BREAK_LAYER_Y_AXIS = true;   % e 图自动省略不同拓扑曲线之间的空白纵轴区间
SHOW_HUB_PANEL = true;
HUB_DEGREE = 6;              % f 图固定目标度，K 与 d 图相同
HUB_N_VALUES = [200];     % f 图只显示 N=200
SHOW_HUB_INSET = true;       % 放大 f<=0.15 的区域
OUTPUT_DIR = '';             % 空值：fig
CACHE_FILE = '';             % 空值：results 中的独立缓存

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir,fullfile(script_dir,'simulations'));
options = struct('code_dir',script_dir,'alpha',ALPHA, ...
    'ratio',COUPLING_RATIO, ...
    'compute',RUN_COMPUTATION,'visible',FIGURE_VISIBLE, ...
    'output_dir',OUTPUT_DIR,'cache_file',CACHE_FILE);

%% 2. 检查已有网络与谱缓存
assert(options.alpha>0 && isfinite(options.alpha) && ...
    options.ratio>=0 && isfinite(options.ratio),'alpha、ratio 无效。');
addpath(options.code_dir);
helper = fullfile(options.code_dir,'multiplex_critical_spectrum.m');
assert(isfile(helper),'缺少已有谱计算函数：%s',helper);
if isempty(options.output_dir)
    options.output_dir=fullfile(options.code_dir,'fig');
end
tag='critical_spectrum_degree_layers';
if isempty(options.cache_file)
    options.cache_file=fullfile(options.code_dir,'results',[tag '.mat']);
end
definitions=makeCases(fullfile(options.code_dir,'results','topology'), ...
    N_VALUES,DEGREE_VALUES,CONNECTIVITY_K,LAYER_VALUES,LAYER_N, ...
    LAYER_DEGREE,INTERLAYER_TYPES,WS_REWIRE);
% 检查已有网络；缓存使用固定文件名，直接复用已有谱点。
files=unique({definitions.topology_file});
available_files=cellfun(@isfile,files);
missing_files=files(~available_files);
for k=1:numel(missing_files)
    fprintf('[跳过] 网络尚未下载：%s\n',missing_files{k});
end
available=cellfun(@isfile,{definitions.topology_file})';
files=files(available_files);
for k=1:numel(files)
    assert(isfile(files{k}),'缺少已有网络；不会重新生成：%s',files{k});
    data=load(files{k},'nets');
    selected=definitions(strcmp({definitions.topology_file},files{k}));
    assert(isfield(data,'nets') && iscell(data.nets) && ...
        numel(data.nets)>=max([selected.K]),'网络文件层数不足：%s',files{k});
end
config=struct('alpha',options.alpha,'ratio',options.ratio, ...
    'schema_version',1);
template=struct('x',NaN,'mean_degree',NaN,'lambda_star',NaN, ...
    'lambda_opt',NaN,'sigma_s',NaN,'growth_at_threshold',NaN, ...
    'method','','status','','residual',NaN,'dense_verified',false);
if isfile(options.cache_file)
    data=load(options.cache_file,'state');
    previous=data.state;
    assert(isequaln(previous.config,config),'谱缓存参数不同，请指定独立 CACHE_FILE。');
    % 按网络和扫描配置匹配缓存，避免选点变化后沿用旧数组下标。
    state=struct('config',config,'definitions',definitions, ...
        'values',repmat(template,numel(definitions),1), ...
        'done',false(numel(definitions),1));
    for idx=1:numel(definitions)
        def=definitions(idx);
        old=find(strcmp({previous.definitions.study},def.study) & ...
            strcmp({previous.definitions.topology},def.topology) & ...
            [previous.definitions.target_degree]==def.target_degree & ...
            [previous.definitions.N]==def.N & [previous.definitions.K]==def.K & ...
            strcmp({previous.definitions.interlayer},def.interlayer),1);
        if ~isempty(old) && previous.done(old)
            state.values(idx)=previous.values(old);
            state.done(idx)=true;
            if strcmp(def.study,'connectivity')
                state.values(idx).x=state.values(idx).mean_degree;
            end
        end
    end
    if isequaln(previous.definitions,definitions) && ...
            isfield(previous,'hub_profiles_by_topology')
        state.hub_profiles_by_topology=previous.hub_profiles_by_topology;
    end
else
    state=struct('config',config,'definitions',definitions, ...
        'values',repmat(template,numel(definitions),1), ...
        'done',false(numel(definitions),1));
end
missing=find(~state.done & available);
assert(options.compute || isempty(missing), ...
    '只读缓存模式缺少 %d 个谱点。',numel(missing));
for idx=reshape(missing,1,[])
    def=definitions(idx);
    fprintf('[谱 %d/%d] %s %s N=%d K=%d\n',idx,numel(definitions), ...
        def.study,def.topology,def.N,def.K);
    data=load(def.topology_file,'nets');
    layers=data.nets(1:def.K);
    degree=zeros(def.K,1);
    for k=1:def.K
        assert(isequal(size(layers{k}),[def.N def.N]),'层内网络尺寸错误。');
        degree(k)=full(mean(diag(layers{k})));
    end
    inter=interlayerLaplacian(def.K,def.interlayer);
    % 低度网络的谱可能整体低于 lambda_opt，允许完整谱回退。
    spectral_options=struct('dense_verify_limit',def.N*def.K);
    spec=multiplex_critical_spectrum(layers,inter,options.alpha,options.ratio, ...
        spectral_options);
    assert(strcmp(spec.status,'ok') && isfinite(spec.lambda_star) && ...
        isfinite(spec.sigma_c_spec),'没有可失稳的非均匀模态：%s',def.topology_file);
    % 模态矩阵在理论临界点应具有零增长率。
    J=[10/3 -5;10 -4];
    modal=J-options.alpha*spec.lambda_star*diag([1 spec.sigma_c_spec]);
    growth=max(real(eig(modal)));
    assert(abs(growth)<1e-7,'模态零增长校验失败：%.3g',growth);
    value=template;
    value.mean_degree=mean(degree);
    switch def.study
        case 'connectivity', value.x=value.mean_degree;
        case 'layer_count', value.x=def.K;
    end
    value.lambda_star=spec.lambda_star;
    value.lambda_opt=spec.lambda_opt;
    value.sigma_s=spec.sigma_c_spec;
    value.growth_at_threshold=growth;
    value.method=spec.method;
    value.status=spec.status;
    value.residual=spec.relative_residual;
    value.dense_verified=spec.dense_verified;
    state.values(idx)=value;
    state.done(idx)=true;
    checkpoint(options.cache_file,state);
end
plotted=state.done & available;
fprintf('[完成] 绘制 %d 个谱点；跳过 %d 个未下载网络对应的谱点。\n', ...
    nnz(plotted),nnz(~available));
if SHOW_HUB_PANEL
    hub_families={'ER','WS','BA'};
    if ~isfield(state,'hub_profiles_by_topology')
        state.hub_profiles_by_topology=struct('ER',struct(),'WS',struct(),'BA',struct());
        % 保留旧版已算好的 BA 权重，按 N 命名，便于以后删减曲线。
        if isfield(state,'hub_profiles')
            cached_definitions=state.definitions(strcmp({state.definitions.study},'connectivity'));
            cached_N=unique([cached_definitions.N],'stable');
            for group=1:min(numel(cached_N),numel(state.hub_profiles))
                if ~isempty(state.hub_profiles{group})
                    key=sprintf('N%d',cached_N(group));
                    state.hub_profiles_by_topology.BA.(key)=state.hub_profiles{group};
                end
            end
        end
    end
    for family=1:numel(hub_families)
    topology=hub_families{family};
    for group=1:numel(HUB_N_VALUES)
        idx=find(plotted' & strcmp({definitions.study},'layer_count') & ...
            strcmp({definitions.topology},topology) & [definitions.N]==HUB_N_VALUES(group) & ...
            [definitions.K]==CONNECTIVITY_K & strcmp({definitions.interlayer},'complete_graph') & ...
            abs([definitions.target_degree]-HUB_DEGREE)<1e-8,1);
        if isempty(idx) && strcmp(topology,'ER')
            candidates=find(plotted' & strcmp({definitions.study},'layer_count') & ...
                strcmp({definitions.topology},topology) & [definitions.N]==HUB_N_VALUES(group) & ...
                [definitions.K]==CONNECTIVITY_K & strcmp({definitions.interlayer},'complete_graph'));
            if ~isempty(candidates)
                [~,nearest]=min(abs([definitions(candidates).target_degree]-HUB_DEGREE));
                idx=candidates(nearest);
            end
        end
        if isempty(idx), continue; end
        key=sprintf('N%d',HUB_N_VALUES(group));
        if isfield(state.hub_profiles_by_topology.(topology),key)
            continue;
        end
        if ~RUN_HUB_PROFILE_COMPUTATION
            fprintf('[跳过] %s N=%d 的临界模式权重尚未缓存。\n',topology,HUB_N_VALUES(group));
            continue;
        end
        def=definitions(idx); data=load(def.topology_file,'nets');
        fprintf('[枢纽] %s N=%d K=%d 目标度=%d\n',topology,def.N,def.K,HUB_DEGREE);
        state.hub_profiles_by_topology.(topology).(key)=critical_mode_degree_profile(data.nets(1:def.K), ...
            interlayerLaplacian(def.K,def.interlayer),options.ratio, ...
            state.values(idx).lambda_star);
        checkpoint(options.cache_file,state);
    end
    end
end
if options.compute,checkpoint(options.cache_file,state);end
if ~isfolder(options.output_dir),mkdir(options.output_dir);end
threshold_config=struct('alpha',ALPHA,'ratio',COUPLING_RATIO, ...
    'noise',THRESHOLD_NOISE,'noise_seed',THRESHOLD_NOISE_SEED, ...
    'T_END',500,'dt',0.001,'steps',2,'init_perturb',0.1, ...
    'detect_convergence',false,'amplitude_cutoff',0.05, ...
    'sigma_tolerance',0.05,'sigma_ceiling',320,'schema_version',1);
threshold_rows=collectTransitionThresholds(state,plotted,script_dir, ...
    threshold_config,RUN_THRESHOLD_COMPUTATION,THRESHOLD_TASK_INDICES,THRESHOLD_WORKERS);
% 与 fig2.m 的 2x3 正式图使用相同窗口尺寸与紧凑布局。
fig=figure('Color','w','Visible',options.visible,'Units','normalized', ...
    'Position',[0.1 0.1 0.75 0.65]);
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
studies={'connectivity','layer_count'};
panel_labels={'(d)','(e)'};
xlabels={'$\langle d\rangle$','$K$'};
families={'ER','WS','BA'};
colors=[0.00 0.447 0.741;0.85 0.325 0.098;0.494 0.184 0.556];
markers={'o','s','^'};
line_styles={'-','--',':'};
interlayer_labels={'R','A'};
[b_axes,b_legend]=plotTransitionPanels(layout,threshold_rows,N_VALUES, ...
    INTERLAYER_TYPES,colors,markers,line_styles,interlayer_labels);
for panel=1:2
    ax=nexttile(layout,panel+3);hold(ax,'on');
    for family=1:3
        if panel==1, groups=numel(N_VALUES); else, groups=numel(INTERLAYER_TYPES); end
        for group=1:groups
        mask=strcmp({definitions.study},studies{panel}) & ...
            strcmp({definitions.topology},families{family}) & plotted';
        if panel==1
            mask=mask & [definitions.N]==N_VALUES(group);
            label=sprintf('%s, N=%d',families{family},N_VALUES(group));
        else
            mask=mask & strcmp({definitions.interlayer},INTERLAYER_TYPES{group});
            label=sprintf('%s-%s',families{family},interlayer_labels{group});
        end
        selected=find(mask);
        if isempty(selected), continue; end
        [x,order]=sort([state.values(selected).x]);
        y=[state.values(selected).lambda_star];y=y(order);
        plot(ax,x,y,'LineStyle',line_styles{group},'Marker',markers{family}, ...
            'Color',colors(family,:),'MarkerFaceColor',colors(family,:), ...
            'MarkerEdgeColor',colors(family,:),'LineWidth',1.8, ...
            'MarkerSize',6,'DisplayName',label);
        end
    end
    xlabel(ax,xlabels{panel},'Interpreter','latex','FontSize',18,'FontName','Times New Roman');
    ylabel(ax,'$\lambda_*$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
    if ~isempty(findobj(ax,'Type','line'))
        if strcmp(panel_labels{panel},'(e)')
            e_legend=legend(ax,'Location','northwest','FontSize',14,'Interpreter','latex','FontName','Times New Roman','NumColumns',3);
            e_axes=ax;
        elseif strcmp(panel_labels{panel},'(d)')
            legend(ax,'Location','southeast','FontSize',14,'Interpreter','latex','FontName','Times New Roman','NumColumns',2);
        else
            legend(ax,'Location','northeast','FontSize',10,'Interpreter','latex','FontName','Times New Roman','NumColumns',2);
        end
    end
    grid(ax,'on');box(ax,'on');
    set(ax,'FontName','Times New Roman','FontSize',18,'TickLabelInterpreter','latex');
    if panel==1
        xlim(ax,[1.8 20]);xticks(ax,2:2:20);
    else
        xlim(ax,[3 10]);xticks(ax,3:10);
        if BREAK_LAYER_Y_AXIS, breakLayerYAxis(ax,families); end
    end
    placePanelLabel(ax,panel_labels{panel},-0.16);
end
if SHOW_HUB_PANEL
ax=nexttile(layout,6);hold(ax,'on');box(ax,'on');xlim(ax,[0 1]);ylim(ax,[0 1]);
xlabel(ax,'$f$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
ylabel(ax,'$C(f)$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
placePanelLabel(ax,'(f)');
curve_count=0;
hub_curves=struct('f',{},'cumulative',{},'color',{},'style',{});
for family=1:numel(families)
topology=families{family};
for group=1:numel(HUB_N_VALUES)
    idx=find(plotted' & strcmp({definitions.study},'layer_count') & ...
        strcmp({definitions.topology},topology) & [definitions.N]==HUB_N_VALUES(group) & ...
        [definitions.K]==CONNECTIVITY_K & strcmp({definitions.interlayer},'complete_graph') & ...
        abs([definitions.target_degree]-HUB_DEGREE)<1e-8,1);
    if isempty(idx) && strcmp(topology,'ER')
        candidates=find(plotted' & strcmp({definitions.study},'layer_count') & ...
            strcmp({definitions.topology},topology) & [definitions.N]==HUB_N_VALUES(group) & ...
            [definitions.K]==CONNECTIVITY_K & strcmp({definitions.interlayer},'complete_graph'));
        if ~isempty(candidates)
            [~,nearest]=min(abs([definitions(candidates).target_degree]-HUB_DEGREE));
            idx=candidates(nearest);
        end
    end
    key=sprintf('N%d',HUB_N_VALUES(group));
    if isempty(idx) || ~isfield(state.hub_profiles_by_topology.(topology),key)
        continue;
    end
    profile=state.hub_profiles_by_topology.(topology).(key);
    if numel(HUB_N_VALUES)==1
        style='-';
    else
        style=line_styles{mod(group-1,numel(line_styles))+1};
    end
    plot(ax,profile.f,profile.cumulative,'Color',colors(family,:), ...
        'LineStyle',style,'LineWidth',1.8, ...
        'DisplayName',sprintf('%s, N=%d',topology,HUB_N_VALUES(group)));
    curve_count=curve_count+1;
    hub_curves(curve_count)=struct('f',profile.f,'cumulative',profile.cumulative, ...
        'color',colors(family,:),'style',style);
    fprintf('[枢纽] %s N=%d：前10%%节点占临界模式权重 %.2f%%，特征空间维数 %d。\n', ...
        topology,HUB_N_VALUES(group),100*profile.top10_weight,profile.multiplicity);
end
end
plot(ax,[0 1],[0 1],'--','Color',[.4 .4 .4], ...
    'LineWidth',1.5,'HandleVisibility','off');
if curve_count==0
    text(ax,.5,.6,'No available critical-mode data','Units','normalized', ...
        'HorizontalAlignment','center','FontSize',21,'Interpreter','latex', ...
        'FontName','Times New Roman','Color',[.45 .45 .45]);
end
legend(ax,'Location','southeast','FontSize',12,'Interpreter','latex','FontName','Times New Roman');grid(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',18,'TickLabelInterpreter','latex');
if SHOW_HUB_INSET && curve_count>0
    drawnow;
    ax.Units='normalized'; position=ax.Position;
    % 放大窗向上移，右下角留给图例。
    zoom_ax=axes(fig,'Units','normalized','Position', ...
        [position(1)+.55*position(3),position(2)+.40*position(4), ...
        .40*position(3),.34*position(4)]);
    hold(zoom_ax,'on');box(zoom_ax,'on');grid(zoom_ax,'on');
    for curve=1:numel(hub_curves)
        curve_data=hub_curves(curve);
        plot(zoom_ax,curve_data.f,curve_data.cumulative,'Color',curve_data.color, ...
            'LineStyle',curve_data.style,'LineWidth',1.5);
    end
    plot(zoom_ax,[0 .15],[0 .15],'--','Color',[.4 .4 .4]);
    xlim(zoom_ax,[0 .15]);ylim(zoom_ax,[0 1]);xticks(zoom_ax,[0 .05 .10 .15]);
    set(zoom_ax,'FontName','Times New Roman','FontSize',10,'TickLabelInterpreter','latex');
end
end

% 最终布局完成后，固定 b、e 图图例在曲线之间的空白区。
if isgraphics(b_legend)
    drawnow;
    b_axes.Units='normalized'; b_legend.Units='normalized';
    axes_position=b_axes.Position; legend_position=b_legend.Position;
    legend_position(1)=axes_position(1)+(axes_position(3)-legend_position(3))/2;
    legend_position(2)=axes_position(2)+0.75*axes_position(4)-legend_position(4);
    b_legend.Position=legend_position;
end
if exist('e_legend','var') && isgraphics(e_legend)
    drawnow;
    e_axes.Units='normalized'; e_legend.Units='normalized';
    axes_position=e_axes.Position; legend_position=e_legend.Position;
    legend_position(1)=axes_position(1)+0.025*axes_position(3);
    legend_position(2)=axes_position(2)+0.72*axes_position(4)-legend_position(4);
    e_legend.Position=legend_position;
end

% 子图标号在 y label 顶部；允许位于坐标框外，但不越过子图外边界。
drawnow;
letters=findall(fig,'Type','text','Tag','panel_letter');
for idx=1:numel(letters)
    letter=letters(idx);ax=letter.Parent;
    ax.Units='normalized';ax.YLabel.Units='normalized';
    label_position=ax.YLabel.Position;
    letter.Position=[label_position(1)-0.025 0.955 0];
end
drawnow;
for idx=1:numel(letters)
    letter=letters(idx);ax=letter.Parent;
    inner=ax.Position;outer=ax.OuterPosition;
    bounds=[(outer(1:2)-inner(1:2))./inner(3:4),outer(3:4)./inner(3:4)];
    extent=letter.Extent;position=letter.Position;
    for dimension=1:2
        low=bounds(dimension);high=low+bounds(dimension+2);
        shift=max(0,low-extent(dimension))-max(0,extent(dimension)+extent(dimension+2)-high);
        position(dimension)=position(dimension)+shift;
    end
    letter.Position=position;
end

% 导出图像：EPS 保存到手稿目录，PNG 保存到本地 fig 目录
test_plots_dir = fullfile(options.code_dir, 'fig');
if ~isfolder(test_plots_dir), mkdir(test_plots_dir); end
png_file = fullfile(test_plots_dir, 'fig5.png');
exportgraphics(fig, png_file, 'Resolution', 300);

base = fullfile(options.output_dir, 'fig5');
exportgraphics(fig, [base '.eps'], 'ContentType', 'vector');

if strcmp(options.visible,'off'),close(fig);end

%% 局部辅助函数
function definitions=makeCases(root,N_values,degrees,connectivity_K, ...
    layer_values,layer_N,layer_degree,interlayer_types,ws_rewire)
definitions=struct('study',{},'topology',{},'N',{},'K',{}, ...
    'target_degree',{},'topology_file',{},'interlayer',{});
families={'ER','WS','BA'};
for N=N_values
    % ER 按已有 p 网格选点，不再由整数目标度反推不存在的文件名。
    er_files=dir(fullfile(root,'ER',sprintf('N%d_p*.mat',N)));
    for f=1:numel(er_files)
        token=regexp(er_files(f).name,'_p([0-9.]+)\.mat$','tokens','once');
        if isempty(token), continue; end
        file=fullfile(er_files(f).folder,er_files(f).name);
        nominal_degree=(N-1)*str2double(token{1});
        % 先缩小候选范围，再按实际平均度筛选。
        if nominal_degree<min(degrees)-2 || nominal_degree>max(degrees)+2
            continue;
        end
        data=load(file,'nets');
        assert(isfield(data,'nets') && iscell(data.nets) && ...
            numel(data.nets)>=connectivity_K,'网络文件层数不足：%s',file);
        actual_degree=mean(cellfun(@(L) full(mean(diag(L))),data.nets(1:connectivity_K)));
        if actual_degree<min(degrees) || actual_degree>max(degrees), continue; end
        definitions(end+1)=item('connectivity','ER',N,connectivity_K, ...
            nominal_degree,file,'complete_graph'); %#ok<AGROW>
    end
    for degree=degrees
        for t=2:3
            [file,supported]=topologyFile(root,families{t},N,degree,ws_rewire);
            if ~supported, continue; end
            definitions(end+1)=item('connectivity',families{t},N, ...
                connectivity_K,degree,file,'complete_graph'); %#ok<AGROW>
        end
    end
end
for K=layer_values
    for t=1:3
        [file,supported]=topologyFile(root,families{t},layer_N,layer_degree,ws_rewire);
        assert(supported,'层数扫描的层内度不受标准网络参数支持。');
        for coupling=1:numel(interlayer_types)
            definitions(end+1)=item('layer_count',families{t},layer_N,K, ...
                layer_degree,file,interlayer_types{coupling}); %#ok<AGROW>
        end
    end
end
end

function [file,supported]=topologyFile(root,family,N,degree,ws_rewire)
% 只读取已有标准库格式；未下载的文件由调用处列出。
supported=true;
switch family
    case 'ER'
        name=sprintf('N%d_p%.3f.mat',N,degree/(N-1));
    case 'WS'
        supported=mod(degree,2)==0;
        name=sprintf('N%d_K%d_pr%.2f.mat',N,degree,ws_rewire);
    case 'BA'
        supported=mod(degree,2)==0;
        name=sprintf('N%d_m%d.mat',N,degree/2);
end
file=fullfile(root,family,name);
end

function def=item(study,topology,N,K,degree,file,interlayer)
def=struct('study',study,'topology',topology,'N',N,'K',K, ...
    'target_degree',degree,'topology_file',file,'interlayer',interlayer);
end

function L=interlayerLaplacian(K,type)
if strcmp(type,'complete_graph')
    L=K*speye(K)-sparse(ones(K));
else
    A=sparse(1:K,[2:K 1],1,K,K);A=A+A';L=2*speye(K)-A;
    if K==2
        fprintf('[层间] K=2 使用权重2连接，以维持加权度2。\n');
    end
end
end

function checkpoint(file,state)
parent=fileparts(file);if ~isfolder(parent),mkdir(parent);end
temporary=[tempname(parent) '.mat'];
save(temporary,'state','-v7');
[ok,message]=movefile(temporary,file,'f');assert(ok,'%s',message);
end


function breakLayerYAxis(ax,families)
% 用分段坐标压缩空白区间；刻度仍显示原始 lambda_* 数值。
curves=findobj(ax,'Type','line');
bands=[];
for family=1:numel(families)
    y=[];
    for j=1:numel(curves)
        if startsWith(curves(j).DisplayName,[families{family} ',']) || ...
                startsWith(curves(j).DisplayName,[families{family} '-'])
            y=[y curves(j).YData(:)']; %#ok<AGROW>
        end
    end
    y=y(isfinite(y));
    if ~isempty(y), bands(end+1,:)=[min(y) max(y)]; end %#ok<AGROW>
end
references=findobj(ax,'Type','constantline');
for j=1:numel(references)
    bands(end+1,:)=[references(j).Value references(j).Value]; %#ok<AGROW>
end
if size(bands,1)<2, return; end
span=max(bands(:,2))-min(bands(:,1));
padding=max(0.15,0.025*span);
bands=sortrows(bands+[-padding padding],1);
merged=bands(1,:);
for j=2:size(bands,1)
    if bands(j,1)<=merged(end,2)+padding
        merged(end,2)=max(merged(end,2),bands(j,2));
    else
        merged(end+1,:)=bands(j,:); %#ok<AGROW>
    end
end
if size(merged,1)<2, return; end
% 每段保留同一比例，仅压缩段与段之间的空白。
gap=padding;
offsets=zeros(size(merged,1),1);
for j=2:size(merged,1)
    offsets(j)=offsets(j-1)+merged(j,1)-merged(j-1,2)-gap;
end
for j=1:numel(curves)
    curves(j).YData=mapBrokenY(curves(j).YData,merged,offsets);
end
for j=1:numel(references)
    references(j).Value=mapBrokenY(references(j).Value,merged,offsets);
end
ticks=[]; labels={};
for j=1:size(merged,1)
    original=linspace(merged(j,1)+padding/3,merged(j,2)-padding/3,3);
    ticks=[ticks original-offsets(j)]; %#ok<AGROW>
    labels=[labels arrayfun(@(v) sprintf('%.2f',v),original,'UniformOutput',false)]; %#ok<AGROW>
end
ylim(ax,[merged(1,1) merged(end,2)-offsets(end)]);
yticks(ax,ticks);yticklabels(ax,labels);
limits=xlim(ax); dx=0.012*diff(limits); dy=0.16*gap;
for j=1:size(merged,1)-1
    low=merged(j,2)-offsets(j); high=merged(j+1,1)-offsets(j+1);
    % 白色横带明确表示省略区间，两侧双斜线标记断轴。
    patch(ax,[limits(1) limits(2) limits(2) limits(1)], ...
        [low low high high],ax.Color,'EdgeColor','none','HandleVisibility','off');
    middle=(low+high)/2;
    for edge=limits
        for shift=[-0.22 0.22]*gap
            line(ax,edge+[-dx dx],middle+shift+[-dy dy], ...
                'Color','k','LineWidth',1,'Clipping','off','HandleVisibility','off');
        end
    end
end
fprintf('[断轴] b 图保留 %d 个纵轴区间，刻度显示原始谱值。\n',size(merged,1));
end

function mapped=mapBrokenY(y,bands,offsets)
mapped=y;
for j=1:size(bands,1)
    selected=y>=bands(j,1) & y<=bands(j,2);
    mapped(selected)=y(selected)-offsets(j);
end
end

function rows=collectTransitionThresholds(source,plotted,script_dir,config,RUN_COMPUTATION,TASK_INDICES,REQUESTED_WORKERS)
definitions=source.definitions;
cache_dir=fullfile(script_dir,'results','topology_transition_thresholds');
if RUN_COMPUTATION && ~isfolder(cache_dir),mkdir(cache_dir);end
cache_index=read_topology_transition_cache_index(cache_dir);
selected=find(plotted);
if ~isempty(TASK_INDICES),selected=intersect(selected,TASK_INDICES,'stable');end
rows=struct('definition',{},'mean_degree',{},'sigma_f',{},'sigma_b',{}, ...
    'forward_bracket',{},'backward_bracket',{},'cache_file',{});
jobs={};
for idx=reshape(selected,1,[])
    def=definitions(idx);
    assert(isfile(def.topology_file),'缺少网络：%s',def.topology_file);
    metadata=struct('definition',def,'config',config);
    [task,cache_file]=load_topology_transition_cache(cache_dir,def,config,cache_index);
    if ~isempty(task)
        task.metadata=metadata;
    else
        task=struct('metadata',metadata,'done',false,'seed',[], ...
            'tests',struct('branch',{},'sigma',{},'amplitude',{}), ...
            'sigma_f',NaN,'sigma_b',NaN,'forward_bracket',[], ...
            'backward_bracket',[]);
    end
    fprintf('[任务 %d/%d] %s %s N=%d K=%d\n',idx,numel(definitions), ...
        def.study,def.topology,def.N,def.K);
    if ~task.done
        if ~RUN_COMPUTATION
            fprintf('[待补算] %s %s N=%d K=%d；仅绘图模式跳过。\n', ...
                def.study,def.topology,def.N,def.K);
            continue;
        end

    end
    jobs{end+1}=struct('definition',def,'value',source.values(idx), ...
        'task',task,'cache_file',cache_file); %#ok<AGROW>
end
results=cell(size(jobs));
pending=cellfun(@(job) ~job.task.done,jobs);
if any(pending)
    validateattributes(REQUESTED_WORKERS,{'numeric'},{'scalar','integer','positive'});
    pending_count=nnz(pending);
    pool=gcp('nocreate');
    if isempty(pool),pool=parpool('Processes',min(REQUESTED_WORKERS,pending_count));end
    workers=min([REQUESTED_WORKERS,pool.NumWorkers,pending_count]);
    fprintf('[并行] %d 个待补算配置，worker=%d。\n',pending_count,workers);
    parfor (job_idx=1:numel(jobs),workers)
        results{job_idx}=computeTransitionJob(jobs{job_idx},config);
    end
else
    for job_idx=1:numel(jobs)
        results{job_idx}=computeTransitionJob(jobs{job_idx},config);
    end
end
if ~isempty(results),rows=[results{:}];end
fprintf('[阈值] 已读取/完成 %d 个配置。\n',numel(rows));
if RUN_COMPUTATION
    summary_file=fullfile(cache_dir,'summary.mat');
    if isfile(summary_file)
        summary_file=[tempname(cache_dir) '_summary.mat'];
    end
    save(summary_file,'rows','config','selected');
end
end

function row=computeTransitionJob(job,config)
def=job.definition;value=job.value;task=job.task;cache_file=job.cache_file;
if ~task.done
        data=load(def.topology_file,'nets');
        assert(iscell(data.nets) && numel(data.nets)>=def.K,'网络层数不足。');
        cfg=config; cfg.N=def.N; cfg.K=def.K; cfg.beta=config.alpha*config.ratio;
        cfg.L_intra=data.nets(1:def.K);
        if strcmp(def.interlayer,'complete_graph')
            cfg.L_inter=def.K*speye(def.K)-sparse(ones(def.K));
        else
            A=sparse(1:def.K,[2:def.K 1],1,def.K,def.K); A=A+A';
            cfg.L_inter=2*speye(def.K)-A;
        end
        for k=1:def.K,assert(isequal(size(cfg.L_intra{k}),[def.N def.N]),'网络尺寸错误。');end
        high=min(config.sigma_ceiling,max(40,1.5*value.sigma_s));
        [low_A,task]=probe(task,cfg,'forward',0,cache_file);
        assert(low_A<=config.amplitude_cutoff,'sigma=0 已形成斑图，无法建立正向阈值区间。');
        while true
            [high_A,task]=probe(task,cfg,'forward',high,cache_file);
            if high_A>config.amplitude_cutoff,break;end
            assert(high<config.sigma_ceiling,'达到 sigma 上限仍无斑图，请检查范围和积分时间。');
            high=min(2*high,config.sigma_ceiling);
        end
        [sigma,bracket,task]=bisect(task,cfg,'forward',high,cache_file);
        task.sigma_f=sigma;task.forward_bracket=bracket;
        % 在同一目标网络上获得强斑图初值，独立保存；不用其他拓扑的状态代替。
        if isempty(task.seed)
            seed_cfg=cfg;seed_cfg.sigma=high;seed_cfg.y0=[];rng(0,'twister');
            [~,Y,out]=solve_multiplex(seed_cfg);
            assert(all(isfinite(Y(end,:))) && isfinite(out.A_final) && ...
                out.A_final>config.amplitude_cutoff,'强斑图初值未形成。');
            task.seed=Y(end,:)';saveTask(cache_file,task);
        end
        [low_A,task]=probe(task,cfg,'backward',0,cache_file);
        [high_A,task]=probe(task,cfg,'backward',high,cache_file);
        assert(low_A<=config.amplitude_cutoff && high_A>config.amplitude_cutoff, ...
            '反向端点未夹住转变，不能报告有效阈值。');
        [sigma,bracket,task]=bisect(task,cfg,'backward',high,cache_file);
        task.sigma_b=sigma;task.backward_bracket=bracket;
        % 二分法要求判据随 sigma 单调；已探测点若相互矛盾则拒绝输出。
        for branch={'forward','backward'}
            tests=task.tests(strcmp({task.tests.branch},branch{1}));
            [~,order]=sort([tests.sigma]);
            patterned=[tests(order).amplitude]>config.amplitude_cutoff;
            assert(all(diff(double(patterned))>=0), ...
                '探测结果非单调，不能用单一二分阈值概括；请检查时间、噪声重复及扫描区间。');
        end
        task.done=true;saveTask(cache_file,task);
end
    row=struct('definition',def,'mean_degree',value.mean_degree, ...
        'sigma_f',task.sigma_f,'sigma_b',task.sigma_b, ...
        'forward_bracket',task.forward_bracket,'backward_bracket',task.backward_bracket, ...
        'cache_file',cache_file);
end

function [b_axes,b_legend]=plotTransitionPanels(layout,rows,N_values, ...
    interlayer_types,colors,markers,line_styles,interlayer_labels)
% 与 d、e 共用样式和组别顺序；缺失某组也不改变其余组的线型。
b_axes=gobjects(0);b_legend=gobjects(0);
panel_labels={'(a)','(b)','(c)'};
families={'ER','WS','BA'};
studies={'connectivity','layer_count','connectivity'};
for panel=1:3
    ax=nexttile(layout,panel);hold(ax,'on');
    for family=1:3
        mask=arrayfun(@(r) strcmp(r.definition.study,studies{panel}) && ...
            strcmp(r.definition.topology,families{family}),rows);
        ids=find(mask);if isempty(ids),continue;end
        if panel~=2
            groups=N_values;
        else
            groups=string(interlayer_types);
        end
        for group=1:numel(groups)
            if panel~=2
                chosen=ids(arrayfun(@(r) r.definition.N==groups(group),rows(ids)));
                x=[rows(chosen).mean_degree]; label=sprintf('%s, N=%d',families{family},groups(group));
            else
                chosen=ids(arrayfun(@(r) string(r.definition.interlayer)==groups(group),rows(ids)));
                x=arrayfun(@(r) r.definition.K,rows(chosen));
                label=sprintf('%s-%s',families{family},interlayer_labels{group});
            end
            if panel==2,branch_separator='-';else,branch_separator=', ';end
            if isempty(chosen),continue;end
            [x,order]=sort(x);chosen=chosen(order);
            if panel==3
                % 直接使用 a 图的配对阈值，不取绝对值、不截断负差值。
                difference=[rows(chosen).sigma_f]-[rows(chosen).sigma_b];
                plot(ax,x,difference,'LineStyle',line_styles{group},'Color',colors(family,:), ...
                    'Marker',markers{family},'MarkerFaceColor',colors(family,:), ...
                    'MarkerEdgeColor',colors(family,:),'LineWidth',1.8,'MarkerSize',6, ...
                    'DisplayName',label);
            else
                % 颜色和标记形状区分网络类型，线型区分组别；实心/空心区分扫描方向。
                plot(ax,x,[rows(chosen).sigma_f],'LineStyle',line_styles{group},'Color',colors(family,:), ...
                    'Marker',markers{family},'MarkerFaceColor',colors(family,:), ...
                    'MarkerEdgeColor',colors(family,:),'LineWidth',1.8,'MarkerSize',6, ...
                    'DisplayName',[label branch_separator 'Fwd']);
                plot(ax,x,[rows(chosen).sigma_b],'LineStyle',line_styles{group},'Color',colors(family,:), ...
                    'Marker',markers{family},'MarkerFaceColor','white', ...
                    'MarkerEdgeColor',colors(family,:),'LineWidth',1.8,'MarkerSize',6, ...
                    'DisplayName',[label branch_separator 'Bwd']);
            end
        end
    end
    if panel~=2,xlabel(ax,'$\langle d\rangle$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
    else,xlabel(ax,'$K$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');end
    if panel==3
        ylabel(ax,'$\Delta\sigma$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
        yline(ax,0,':','Color',[.4 .4 .4],'LineWidth',1.5,'HandleVisibility','off');
    else
        ylabel(ax,'$\sigma$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
    end
    placePanelLabel(ax,sprintf('(%c)','a'+panel-1),-0.16);
    grid(ax,'on');box(ax,'on');
    if isempty(findobj(ax,'Type','line'))
        text(ax,.5,.5,'No cached thresholds','Units','normalized', ...
            'HorizontalAlignment','center','FontSize',21,'Interpreter','latex', ...
            'FontName','Times New Roman','Color',[.5 .5 .5]);
    else
        if strcmp(panel_labels{panel},'(b)')
            b_axes=ax;
            b_legend=legend(ax,'Location','best','FontSize',10,'Interpreter','latex','FontName','Times New Roman','NumColumns',3);
        else
            legend_font_size=10;
            if panel==3,legend_font_size=18;end
            legend(ax,'Location','northeast','FontSize',legend_font_size,'Interpreter','latex','FontName','Times New Roman','NumColumns',2);
        end
    end
    set(ax,'FontName','Times New Roman','FontSize',18,'TickLabelInterpreter','latex');
    if panel==2
        xlim(ax,[3 10]);xticks(ax,3:10);
    else
        xlim(ax,[1.8 20]);xticks(ax,2:2:20);
    end
    if panel==1,ylim(ax,[11 35]);end
end
end

function [sigma,bracket,task]=bisect(task,cfg,branch,high,path)
low=0;
while high-low>cfg.sigma_tolerance
    mid=(low+high)/2;
    [amplitude,task]=probe(task,cfg,branch,mid,path);
    if amplitude>cfg.amplitude_cutoff,high=mid;else,low=mid;end
end
bracket=[low high];sigma=mean(bracket);
end

function [amplitude,task]=probe(task,cfg,branch,sigma,path)
idx=find(strcmp({task.tests.branch},branch) & [task.tests.sigma]==sigma,1);
if ~isempty(idx),amplitude=task.tests(idx).amplitude;return;end
cfg.sigma=sigma;
if strcmp(branch,'forward'),cfg.y0=[];else,cfg.y0=task.seed;end
rng(0,'twister');[~,~,out]=solve_multiplex(cfg);amplitude=out.A_final;
assert(isfinite(amplitude),'积分结果非有限数值，请检查步长与参数。');
task.tests(end+1)=struct('branch',branch,'sigma',sigma,'amplitude',amplitude);
saveTask(path,task);
fprintf('  [%s] sigma=%.5g，A=%.5g\n',branch,sigma,amplitude);
end

function saveTask(path,task)
temporary=[tempname(fileparts(path)) '.mat'];save(temporary,'task','-v7');
[ok,message]=movefile(temporary,path,'f');assert(ok,'%s',message);
end

function placePanelLabel(ax,label,x_position)
% 最终布局完成后，标号与 y label 对齐并限制在所属子图外边界内。
if nargin<3,x_position=-0.12;end
text(ax,x_position,1.01,['\textbf{' label '}'],'Units','normalized', ...
    'FontName','Times New Roman','FontSize',21,'Interpreter','latex', ...
    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'Clipping','off','HandleVisibility','off','Tag','panel_letter');
end
