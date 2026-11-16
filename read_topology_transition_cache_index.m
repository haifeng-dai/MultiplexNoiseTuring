function index=read_topology_transition_cache_index(cache_dir)
% 每次运行只扫描一次；同配置的多个文件保留供读取时检查冲突。
index=containers.Map('KeyType','char','ValueType','any');
files=dir(fullfile(cache_dir,'*.mat'));
loaded=0;
for idx=1:numel(files)
    path=fullfile(files(idx).folder,files(idx).name);
    variables=whos('-file',path);
    if ~any(strcmp({variables.name},'task')),continue;end
    saved=load(path,'task');task=saved.task;
    if ~isfield(task,'metadata') || ~isfield(task.metadata,'definition') || ...
            ~isfield(task.metadata,'config'),continue;end
    key=topology_transition_cache_key(task.metadata.definition,task.metadata.config);
    entry=struct('path',path,'task',task);
    if isKey(index,key)
        entries=index(key);entries{end+1}=entry;index(key)=entries;
    else
        index(key)={entry};
    end
    loaded=loaded+1;
end
fprintf('[缓存索引] 扫描 %d 个 MAT 文件，加载 %d 个任务缓存，%d 个配置。\n', ...
    numel(files),loaded,index.Count);
end
