function key=topology_transition_cache_key(definition,config)
% 使用配置文本作为内存索引，不计算哈希，不包含机器绝对路径。
[~,network]=fileparts(definition.topology_file);
identity=struct('study',definition.study,'topology',definition.topology, ...
    'N',definition.N,'K',definition.K,'target_degree',definition.target_degree, ...
    'interlayer',definition.interlayer,'network',network,'config',orderfields(config));
key=jsonencode(identity);
end
