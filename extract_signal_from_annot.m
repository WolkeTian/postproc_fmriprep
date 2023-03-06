function [label_name, timeseries] = extract_signal_from_annot(giidata, annotfile, method)
    % L 每个顶点对应结构的索引; ct.table(:,5); 每个结构对应的索引
    % method: pca or mean
    
    [~, L, ct]  = read_annotation(annotfile);
    label_index = ct.table(:,5);
    label_name = ct.struct_names;
    timeseries = zeros(size(giidata,2), numel(label_name)); % times * labels

    if (nargin<3)
        method = 'mean';  % 如果未出现该变量，则对其进行赋值
    end
    disp(['信号提取方式为', method]);

    for x = 1:numel(label_name)
        % 获取第x个label对应的顶点的时间序列
        parcel_data = giidata(L == label_index(x), :);
        
        if size(parcel_data,1) 
            % 如果提取出有数据; 部分模板的label_index可能包含不对应L的值
            if isequal(method,'mean')
                % 获取结构信号均值
                parcel_signal = mean(parcel_data)';
            elseif isequal(method,'pca')
                %parcel_signal = pca(parcel_data); % 此为主成分系数
                %parcel_signal = parcel_signal(:,1);
    
                [~,score] = pca(parcel_data');
                parcel_signal = score(:,1);
    
            end
        else
            % 提取出空数据, 无论方法, 用nan值（mean得到）
            parcel_signal = mean(parcel_data)'; 
        end
        timeseries(:, x) = parcel_signal;
        
    end
    
end