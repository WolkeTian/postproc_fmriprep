function noise_vars = confounds_select(confounds_tsv)
        % load covariates
        % 输入变量为dir得到的结构体
        %noiseIC_index = load(fullfile(AROMA_IC_csv.folder, AROMA_IC_csv.name)); % IC AROMA的头动噪声变量
        confounds_all = tdfread(fullfile(confounds_tsv.folder, confounds_tsv.name)); % read tsvdata
        %aroma_data = confounds_all.arom
        fields = fieldnames(confounds_all); % 搜索filedname
        aroma_index = strfind(fields, 'aroma');
        aroma_index = cellfun(@(x) ~isempty(x), aroma_index); % 转为0;1逻辑向量;1代表是aroma变量
        disp(['基于ICA_aroma, 被归为噪声成分的数量为',num2str(sum(aroma_index))]);
        % motion异常值/箱形函数的逻辑变量(0,1索引）
        % motion_outlier = strfind(fields, 'motion_outlier'); motion_outlier = cellfun(@(x) ~isempty(x), motion_outlier);
        %disp(['基于motion异常点, 回归的数量为',num2str(sum(motion_outlier))]);
        % 不稳定状态(0,1索引）
        unsteady = strfind(fields, 'non_steady_state_outlier'); unsteady = cellfun(@(x) ~isempty(x), unsteady); 
        if sum(unsteady) == 0
            % 如果没有unsteady_state_outlier
            unsteady = zeros(size(aroma_index,1), 1);
        end


        % 高通滤波(0,1索引）
        highpass = strfind(fields, 'cosine'); highpass = cellfun(@(x) ~isempty(x), highpass); 

        %% 如果要启用acompcor, 则启用下方列
        %acompcor 50% variance combined mask
        acompcor_index = strfind(fields, 'a_comp_cor');
        %得到列的逻辑索引
        acompcor_index = cellfun(@(x) ~isempty(x), acompcor_index);
        disp(['基于combined mask, 解释50%方差的acompcor成份数为',num2str(sum(acompcor_index))]);

        % 得到c_compcor_index和w_compcor_index
        % c_compcor保留所有导出的成分，即解释50%方差的ccompcor成分
        ccompcor_index = strfind(fields, 'c_comp_cor');
        ccompcor_index = cellfun(@(x) ~isempty(x), ccompcor_index);
        % w_compcor保留所有导出的成分，即解释50%方差的ccompcor成分
        wcompcor_index = strfind(fields, 'w_comp_cor');
        wcompcor_index = cellfun(@(x) ~isempty(x), wcompcor_index);
% 
%         %%%% 如果w_compcor仅保留前5个成分，则启用下行
%         temp_location = find(wcompcor_index);
%         temp_location = temp_location(1:5);
%         wcompcor_index = zeros(size(ccompcor_index));
%         wcompcor_index(temp_location) = 1; % 仅保留前5

% 
%         disp(['acompcor成份数为白质',num2str(sum(wcompcor_index)), '和脑脊液', num2str(sum(ccompcor_index))]);
%%       acompcor相关列结束


%         confounds_compcor = [confounds_all.c_comp_cor_00, confounds_all.c_comp_cor_01, ...
%             confounds_all.c_comp_cor_02, confounds_all.c_comp_cor_03, confounds_all.c_comp_cor_04, ...
%             confounds_all.w_comp_cor_00, confounds_all.w_comp_cor_01, confounds_all.w_comp_cor_02, ...
%             confounds_all.w_comp_cor_03, confounds_all.w_comp_cor_04];
        %获取aroma，acompcor,不稳定状态和motion异常
        %noise_indexs_all = find(motion_outlier + unsteady + aroma_index + acompcor_index + highpass); % 在噪声变量文件中的列

        % 如果是acompcor 50% 
        % disp('按照aroma + acompcor 50% + highpass去除噪声')
        %noise_indexs_all = find(unsteady + aroma_index + acompcor_index + highpass);

        % 如果是ccompcor+wcompcor_index 50% 
        % disp('按照aroma + ccompcor+wcompcor_index 50% + highpass去除噪声')
        % noise_indexs_all = find(unsteady + aroma_index + wcompcor_index + ccompcor_index + highpass);

        % 如果是仅aroma 50%，也不滤波
        disp('近按照ICA-aroma和去除unsteady点')
        noise_indexs_all = find(unsteady + aroma_index);


        noise_all = zeros(size(confounds_all.global_signal, 1), numel(noise_indexs_all));
        for n = 1:numel(noise_indexs_all)
            noise_all(:, n) = eval(['confounds_all.', fields{noise_indexs_all(n)}]);
        end

        %避免unsteady和motion_outlier有完全一样的列--相关为1
        corrmat = tril(corr(noise_all), -1); [a,~] = find(round(corrmat,5) == 1);
        % 去除其中多余重复的列
        noise_all(:,a) = [];
        
        % 组合所有目标噪声变量
        % noises_2regress = [confounds_compcor, noise_all];
        noises_2regress = noise_all;
        % 存储目标噪声变量，供task使用
        matfname = fullfile(confounds_tsv.folder, [confounds_tsv.name(1:end-4), '_noise_selected.mat']);
        save(matfname, 'noises_2regress');
        % 得到output
        disp(['共回归噪声变量列数为',num2str(size(noises_2regress, 2))]);
        noise_vars = noises_2regress;
end
