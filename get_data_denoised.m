function data_denoised = get_data_denoised(data2D, noise_vars)
    % 从时间*维度的2维数据去除协变量
    disp('——————————————数据去噪中——————————————————');
    mask = logical(std(data2D)); % 标准差为0的不纳入计算
    data_2cal = data2D(:, mask);
    % 保留每个体素的中位数信号强度值(排除前10个点，防止开始时可能存在不稳定的超高值)，后续加回来
    data_mean = repmat(mean(data_2cal(11:end,:), 1), size(data2D,1), 1); 
    %%% 替代方案为保留中位数强度值
    %data_mean = repmat(median(data_2cal, 1), size(data2D,1), 1); 
    data_residual = reg_corr(data_2cal, noise_vars);
    data_cleaned = data_residual + data_mean;

    data2D(:, mask) = data_cleaned; % mask内的重新赋值为去噪后的数据
    data_denoised = data2D;
end