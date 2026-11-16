function [task,path]=load_topology_transition_cache(cache_dir,definition,config,index)
% 使用可读文件名；兼容旧缓存，不依赖机器路径或代码/文件哈希。
[~,network]=fileparts(definition.topology_file);
name=sprintf('%s_%s_%s_K%d_%s_a%.12g_r%.12g_noise%.12g_seed%d.mat', ...
    definition.study,definition.topology,network,definition.K, ...
    definition.interlayer,config.alpha,config.ratio,config.noise,config.noise_seed);
path=fullfile(cache_dir,name);
task=[];
key=topology_transition_cache_key(definition,config);
if ~isKey(index,key),return;end
entries=index(key);
for idx=1:numel(entries)
    candidate=entries{idx}.path;
    old=entries{idx}.task;
    % 旧版多输出赋值覆盖了阈值字段；用保存的探测结果恢复二分区间。
    for branch={'forward','backward'}
        if strcmp(branch{1},'forward'),field='sigma_f';else,field='sigma_b';end
        if isfinite(old.(field)),continue;end
        tests=old.tests(strcmp({old.tests.branch},branch{1}));
        amplitudes=[tests.amplitude];sigmas=[tests.sigma];
        valid=isfinite(amplitudes) & isfinite(sigmas);
        low=sigmas(valid & amplitudes<=config.amplitude_cutoff);
        high=sigmas(valid & amplitudes>config.amplitude_cutoff);
        if isempty(low) || isempty(high),continue;end
        bracket=[max(low) min(high)];
        if bracket(2)>=bracket(1) && diff(bracket)<=config.sigma_tolerance*(1+1e-10)
            old.(field)=mean(bracket);
            old.([branch{1} '_bracket'])=bracket;
        end
    end
    old.done=old.done && isfinite(old.sigma_f) && isfinite(old.sigma_b);
    if isempty(task)
        task=old;
    elseif ~task.done && old.done
        task=old;
    elseif task.done && old.done
        assert(isequaln(task.sigma_f,old.sigma_f) && isequaln(task.sigma_b,old.sigma_b), ...
            '同一配置存在不同的已完成阈值，请检查缓存：%s',candidate);
    elseif ~task.done && ~old.done && numel(old.tests)>numel(task.tests)
        task=old;
    end
end
end
