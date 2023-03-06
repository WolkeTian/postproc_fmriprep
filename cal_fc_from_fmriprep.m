function cal_fc_from_fmriprep(bidsout,method)
    % store output text to diary file
    diaryfile = ['Postproc_diary_', char(datetime), '.txt'];
    diaryfile = fullfile(bidsout, strrep(diaryfile, ':', '_'));
    % create diary text file
    fclose(fopen(diaryfile, 'a'));
    diary(diaryfile);
    %% 开始运行

    disp(['extract time series and calculate resting-state FC from ', bidsout]);
    
    % 移除运动异常点，ICA-AROMA对高运动已有不错的效果，不需要censor; 实测增加与否影响很小，提取出的信号平均相关0.98以上
    % acompcor 50 和c/wcompcor5差别还是大点，实测大致相关在0.8（0.6-1之间）;对最后得到功能连接的差异略小些。
    % c/wcompcor 5 和 c/wcompcor 50 差别也差不多，平均相关大概在0.8出头;对于皮层和皮层下数据情况类似。
    % ica-aroma与ica-aroma结合acompcor: 平均0.82左右相关
    %          与结合acompcor50相比，平均0.67左右相关
    % acompcor50加ica-aroma, 与没有ica-aroma相比，平均相关0.7左右且方差较大。
    % 仅ica-aroma与没有ica-aroma相比, 信号平均相关仅0.45左右，说明运动噪声影响确实很大。
    % 仅acompcor50与没有相比, 信号平均相关0.66左右，方差很大(0.5几-0.7几)
    % 仅使用acompcor50和仅使用aroma相比，平均相关0.53左右。
    % unsteady_state对相关计算影响挺大,竟然平均相关只有0.4甚至0.1
    % 包含unsteady_state会对结果有很大的影响; 但不回归则连接值无法反映真实值(很虚高);
    % 回归unsteady_state后和删除时间点计算出来的连接值相同（等效）
    % highpass对FC计算影响较小, 微弱降低方差, 连接矩阵相关0.9以上
    % 不同次/fmriprep版本run的预处理得出来的值也略有差异。
    % 回归的协变量越多，方差越小;损失自由度越高。需要权衡
    % 对比sc-fc相关: 基线(0.1相关左右) < ICA-aroma < aroma+acompcor5 < aroma +
    % w/c_compcor5 < aroma + acompcor50或w/c_compcor50;pca提取出结果不如mean提取
    % (可能提升到0.14左右,越到后边提升越小)
    % 最终决定回归ica-aroma+acompcor50(combined-mask)
    disp('默认回归不稳定状态时间点');
    disp('       Discrete cosine-basis regressors, ICA-AROMA噪声IC');
    disp('       以及acompcor基于白质和脑脊液组织combined-mask（包含50%方差成分数量）');
    disp('去噪参考文献：https://www.sciencedirect.com/science/article/pii/S1053811917310972')
    disp('如需调整请修改代码');
    disp('依赖spm12');
    addpath(genpath('.\NifTI_20140122')); % 临时添加nifti文件夹
    %addpath(genpath('.\atlases_surfaces'));
    addpath(genpath('.\freesurfer'));
    % method: pca or mean
    if (nargin<2)
        method = 'mean';  % 提取roi内信号的方法;如果未出现该变量，则对其进行赋值
    end
    disp(['提取ROI内信号的方法为', method]);
    tic;

    % 搜索协变量文件, 如有则考虑分析此func文件夹下的文件
    func_getconfounds = dir(fullfile(bidsout, 'sub-*', 'ses-*', 'func', 'sub*task-rest*desc-confounds*timeseries*tsv'));
    if numel(func_getconfounds) == 0 % 如果没搜到则去掉ses标志
        func_getconfounds = dir(fullfile(bidsout, 'sub-*', 'func', 'sub*task-rest*desc-confounds*timeseries*tsv'));
    end
    
    subcortical_labels = tdfread(fullfile(bidsout, 'desc-aseg_dseg.tsv')); % 读取皮层下label信息
    disp('开始正式运算');
    for i = 1:numel(func_getconfounds)
        func_folder = func_getconfounds(i).folder; % 于当前文件夹下进行分析
        disp(['当前处理路径为',func_folder]);
        % 创建储存FC的文件夹
        FC_folder = fullfile(func_folder, ['connectome_', method,'_signal']);
        mkdir(FC_folder);
        %% 准备协变量部分
        %__________________________________________________________________
        % covariates
        %AROMA_IC_csv = dir(fullfile(func_folder, 'sub*task-rest*AROMAnoiseICs*csv'));
        confounds_tsv = dir(fullfile(func_folder, 'sub*task-rest*desc-confounds*timeseries*tsv'));
        % get covariates,传入结构体
        % 编辑目标协变量需更改confounds_select函数
        noise_vars = confounds_select(confounds_tsv);
        %% 准备volume文件及去噪
        %__________________________________________________________________
        % prepare preprocessed nii files to extract subcortical data
        preproc_nii = dir(fullfile(func_folder, 'sub*task-rest*MNI152*preproc*bold*.nii.gz'));
        gunzip(fullfile(preproc_nii.folder, preproc_nii.name));
        niiname = fullfile(preproc_nii.folder, preproc_nii.name(1:end-3)); % get nifti fnames
        niftiobj = load_nii(niiname);
        % 删除文件
        delete(niiname);
        niftidata = niftiobj.img;
        niftidata_2d = reshape(niftidata, [size(niftidata,1) * size(niftidata,2) * size(niftidata,3),...
            size(niftidata,4)])'; % 转为2d, 第一维为时间
        data_denoised = get_data_denoised(niftidata_2d, noise_vars);
        disp(['当前处理文件为',preproc_nii.name]);
        data_denoised = data_denoised';
        % 得到去噪后的volume数据
        niftidata_denoised = reshape(data_denoised, size(niftidata));
        % 保留nifti文件
        niftiobj.img = niftidata_denoised;
        new_niftiname = [replace(niftiobj.fileprefix, 'desc-preproc', 'desc-denoised'), '.nii'];
        niftiobj.fileprefix = replace(niftiobj.fileprefix, 'desc-preproc', 'desc-denoised');
        niftiobj.hdr.dime.scl_inter = 0; % 避免图像基础值增加显示异常
        save_nii(niftiobj, new_niftiname);
        gzip(new_niftiname); delete(new_niftiname);
        %% get subcortical data,提取皮层下切割后的时间序列
        %__________________________________________________________________
        % from niftidata_denoised
        submask = dir(fullfile(func_folder, 'sub*task-rest*MNI152*desc-aseg_dseg*nii.gz'));
        gunzip(fullfile(submask.folder, submask.name));
        sub_niifile = fullfile(submask.folder, submask.name); 
        sub_niifile = sub_niifile(1:end-3); % 获取皮层下分割的nii名
        submask_obj = load_nifti(sub_niifile);
        % 删除解压后的皮层下mask文件
        delete(sub_niifile);
        submask_values = submask_obj.vol(:); % 转为一维
        % 转为2d,且排除0/unkonw区域
        submask_index = unique(submask_values); submask_index = submask_index(submask_index ~= 0);
        % 初始化变量，存储下皮层信号和label，供功能连接计算
        subsignal_2use = zeros(size(noise_vars,1), numel(submask_index));
        subsignal_2use_labels = cell(numel(submask_index), 1);
        for m = 1:numel(submask_index)
            % 找到mask文件中roi索引对应的位置
            order_in_dseg = find(subcortical_labels.index == submask_index(m));
            subcortical(m).label = subcortical_labels.name(order_in_dseg,:);
            subsignal_2use_labels{m} = subcortical(m).label;
            % 提取平均时间序列
            subcor_signal = data_denoised(submask_values == submask_index(m), :);
            
            if isequal(method, 'mean')
                % 计算均值
                subcortical(m).timeseries = mean(subcor_signal)';
            elseif isequal(method, 'pca')
                % 提取pca第一主成分
                [~,score] = pca(subcor_signal');
                subcortical(m).timeseries = score(:,1);
                
            end
            subsignal_2use(:, m) = subcortical(m).timeseries;
        end
        % 保存提取出来的subcortical timeseries
        matname = fullfile(FC_folder,'Denoised_timeseries_subcortical.mat');

        save(matname, 'subcortical');
        clear subcortical
        % 计算皮层下相关矩阵并保存
        % 计算相关矩阵
        r_mat = corr(subsignal_2use);
        r_mat = tril(r_mat, -1) + tril(r_mat, -1)'; % 中间变为0
        % fisher-z转换
        FCmat = atanh(r_mat);
        % 保存相关矩阵,仅为皮层间
        connect_matname = fullfile(FC_folder,'Subcortical_connectome.mat');
        save(connect_matname, 'FCmat'); 
        % 保存label
        label_name = replace(connect_matname, '_connectome.mat', '_label.mat');
        save(label_name, 'subsignal_2use_labels');

        %% prepare fsaverage surface gii files 去噪并保存去噪后gii文件
        %__________________________________________________________________
        giifile_LH = dir(fullfile(func_folder, 'sub*task-rest*hemi-L*fsaverage*bold*func.gii'));
        giifile_RH = dir(fullfile(func_folder, 'sub*task-rest*hemi-R*fsaverage*bold*func.gii'));
        % load gii file
        giistruct_LH = gifti(fullfile(giifile_LH.folder, giifile_LH.name));
        giidata_LH = giistruct_LH.cdata; % 163842*205的2维数据
        % 去噪
        giidata_LH_denoised = get_data_denoised(giidata_LH', noise_vars);
        disp(['当前处理文件为',giifile_LH.name]);
        giidata_LH_denoised = giidata_LH_denoised';
        % load gii file
        giistruct_RH = gifti(fullfile(giifile_RH.folder, giifile_RH.name));
        giidata_RH = giistruct_RH.cdata;
        %去噪
        giidata_RH_denoised = get_data_denoised(giidata_RH', noise_vars);
        disp(['当前处理文件为',giifile_RH.name]);
        giidata_RH_denoised = giidata_RH_denoised';
        
        % 保存,使用spm12/@gifti下的函数
        % 替换关键词，防止非第一次运行时giifile_LH结构体检索出多个文件
        new_giiname_LH = replace(fullfile(giifile_LH.folder, giifile_LH.name), '.func', '_desc-denoised');  % 避开func(也是文件夹名）
        new_giiname_RH = replace(fullfile(giifile_RH.folder, giifile_RH.name), '.func', '_desc-denoised'); 

        %giistruct_LH.cdata = giidata_LH_denoised; % 不能直接数据替换
        new_gii_struct = gifti(giidata_LH_denoised);
        save(new_gii_struct, new_giiname_LH);
        % 保存右脑
        new_gii_struct = gifti(giidata_RH_denoised);
        save(new_gii_struct, new_giiname_RH);

        %% 开始读取annot皮层模板数据,并计算连接矩阵
        %__________________________________________________________________

        annots_lh = dir('atlases_surfaces\lh*.annot');
        annots_rh = dir('atlases_surfaces\rh*.annot');
        disp('开始计算功能连接');
        for xxx = 1:numel(annots_lh)
            disp(['当前处理皮层分割文件为 ',annots_lh(xxx).name]);
            annotfile = fullfile('atlases_surfaces', annots_lh(xxx).name);
            [label_LH, timeseries_LH] = extract_signal_from_annot(giidata_LH_denoised, annotfile, method); %按照均值提取
            annotfile = fullfile('atlases_surfaces', annots_rh(xxx).name);
            [label_RH, timeseries_RH] = extract_signal_from_annot(giidata_RH_denoised, annotfile, method); %按照均值提取
            % 保存至变量
            LH_signal(xxx).label = label_LH; LH_signal(xxx).timeseries = timeseries_LH;
            RH_signal(xxx).label = label_RH; RH_signal(xxx).timeseries = timeseries_RH;

            
            % 计算功能连接 待完成
            % 丢弃掉第一个label和时间序列,一般为unknown或背景
            label_LH_2use = label_LH(2:end); signal_LH_2use = timeseries_LH(:, 2:end);
            label_RH_2use = label_RH(2:end); signal_RH_2use = timeseries_RH(:, 2:end);

            surf_label_all = [label_LH_2use; label_RH_2use];
            surf_signal_2use = [signal_LH_2use, signal_RH_2use];
            % 计算相关矩阵
            r_mat = corr(surf_signal_2use);
            r_mat = tril(r_mat, -1) + tril(r_mat, -1)'; % 中间变为0
            % fisher-z转换
            FCmat = atanh(r_mat);
            % 保存相关矩阵,仅为皮层间
            tempname = fullfile(FC_folder,annots_lh(xxx).name);
            tempname = replace(tempname, '\lh.', '\');
            connect_matname = replace(tempname, '.annot', '_surface_connectome.mat');
            save(connect_matname, 'FCmat'); 
            % 保存label
            label_name = replace(connect_matname, '_connectome.mat', '_label.mat');
            save(label_name, 'surf_label_all');

            %% 包含皮层下数据，计算连接矩阵
            full_labels = [surf_label_all; subsignal_2use_labels];
            full_signal = [surf_signal_2use, subsignal_2use];
            % 计算相关矩阵
            r_mat = corr(full_signal);
            r_mat = tril(r_mat, -1) + tril(r_mat, -1)'; % 中间变为0
            % fisher-z转换
            FCmat = atanh(r_mat);
            % 保存相关矩阵,仅为皮层间
            
            connect_matname = replace(tempname, '.annot', '_full_connectome.mat');
            save(connect_matname, 'FCmat'); 
            % 保存label
            label_name = replace(connect_matname, '_connectome.mat', '_label.mat');
            save(label_name, 'full_labels');

        end
        matname = fullfile(FC_folder,'Denoised_timeseries_surf_from_annots.mat');
        %save matname LH_signal RH_signal

        save(matname, 'LH_signal', 'RH_signal');
        clear LH_signal RH_signal; % 清除结构体
        disp([func_folder,'当前处理路径已完成']);
    end
    rmpath(genpath('.\NifTI_20140122')); % 临时添加nifti文件夹
    %rmpath(genpath('.\atlases_surfaces'));
    rmpath(genpath('.\freesurfer'));
    %% 计算全部完成
    % 保存日志
    task.diary = onCleanup(@() diary('off')); % when error to stop diary;
    disp('--------------------------------------------------------');
    disp('---------The process for all has been finished---------');
    disp('--------------------------------------------------------');

    toc;diary off;

end