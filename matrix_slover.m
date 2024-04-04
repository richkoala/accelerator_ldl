classdef matrix_slover
    methods  
        %% 稀疏矩阵Lx = b 求解x 
        function [x] = slove_lx(obj,N,LD,b)

             b_t = b;               %% 待处理的右边项
             
             for j=1:N
                nnz_colj = length(find(LD(:,2)==j));   % 当前列非零元素数量
                idx_p1   = find(LD(:,2)==j);           % 当前列 所有非零元素 对应的LD矩阵中位置（向量）
                
                for i=2:nnz_colj
                    r_idx_m =  LD(idx_p1(i),1);        % 待修正的行序号

                    
                    idx_m   = idx_p1(i);               % 待修正元素在稀疏矩阵LD中的位置            
                    
                    %b_i = b_i - b_j*x_ij
                    b_t(r_idx_m) = b_t(r_idx_m) - b_t(j)*LD(idx_m,3);
                end
                
             end
             
             x = b_t;
            
        end
        
        %% Dx = b 求解x
        function [x] = slove_dx(obj,N,LD,b)

             b_t = b;               %% 待处理的右边项
             
             for j=1:N
                idx_p1   = find(LD(:,2)==j);    % 当前列 所有非零元素 对应的LD矩阵中位置（向量）
                idx_m = idx_p1(1);              % 当前列的第一个元素对应的序号即为D元素在LD矩阵中位置             
                i = j;
                %b_i = b_i / x_jj
                b_t(i) = b_t(i) / LD(idx_m,3);
                
             end
             
             x = b_t;
            
        end
        
        %% 稀疏矩阵Lx = b 求解x 
        function [x] = slove_ltx(obj,N,LD,b)

             b_t = b;               %% 待处理的右边项
             
             for i = N:-1:1         %% L转置后上三角矩阵，需要从下向上求解处理，即需要从LD矩阵的第i行向上求解
                 
                nnz_coli = length(find(LD(:,1)==i));   % 当前列非零元素数量
                idx_p1   = find(LD(:,1)==i);           % 当前列 所有非零元素 对应的LD矩阵中位置（向量）
                
                for j = 1:nnz_coli-1
                    r_idx_m =  LD(idx_p1(j),2);        % 待修正的行序号,需要考虑转置，应该是LD中的列序号

                    
                    idx_m   = idx_p1(j);               % 待修正元素在稀疏矩阵LD中的位置            
                    
                    %b_i = b_i - b_j*x_ij
                    b_t(r_idx_m) = b_t(r_idx_m) - b_t(i)*LD(idx_m,3);
                end
                
             end
             
             x = b_t;
            
        end
        
        
        %% 稀疏矩阵Lx = b 求解x ， 引入优先级 
        function [x,opt] = slove_lx_prior(obj,N,LD,b)
             CMD_SOLVE_LX = 3;
             CMD_SOLVE_LTX= 4;
             CMD_SOLVE_DX = 5;
             
             prior = ones(N,1);     %% 设置所有的优先级均为1，由于L矩阵数据固定，只有b变更存在依赖关系，       
             
             %统计总的运算量
             opt_num = length(LD) - N;
             opt     = zeros(opt_num,6);
             opt_idx = 0;
             
             b_t = b;               %% 待处理的右边项
             
             for j=1:N
                nnz_colj = length(find(LD(:,2)==j));   % 当前列非零元素数量
                idx_p1   = find(LD(:,2)==j);           % 当前列非零元素对应的行序号
                
                for i=2:nnz_colj
                    r_idx_m =  LD(idx_p1(i),1);        % 待修正的行序号                   
                    idx_m   = idx_p1(i);               % 待修正元素在稀疏矩阵LD中的位置            
                    
                    %b_i = b_i - b_j*x_ij
                    b_t(r_idx_m) = b_t(r_idx_m) - b_t(j)*LD(idx_m,3);
                    
                    idx0 = r_idx_m;         % b向量列被更新元素行序号
                    idx1 = j;               % b向量列更新需要使用b向量中元素的行序号
                    idx2 = idx_m;           % b向量列更新需要使用LD矩阵中元素中序号（1维）
                    idx3 = 1;               % 没有使用
                   
                    prior(r_idx_m)   = max(prior(r_idx_m),prior(j))+1;
                    
                    opt_idx          = opt_idx+1;
                    opt(opt_idx,1)   = prior(r_idx_m);
                    opt(opt_idx,2)   = CMD_SOLVE_LX;
                    opt(opt_idx,3:6) = [idx0,idx1,idx2,idx3];
                    
                end
                
             end
             
             x = b_t;
             
             % opt_num - opt_num_th %理论计算与统计的偏差
            
        end
        
        
        function [x,opt] = slove_dx_prior(obj,N,LD,b)
            
             CMD_SOLVE_LX = 3;
             CMD_SOLVE_LTX= 4;
             CMD_SOLVE_DX = 5;

             b_t = b;               %% 待处理的右边项
             opt_num = N;
             prior = ones(N,1);     %% 设置所有的优先级均为1 
             
             opt     = zeros(opt_num,6);
             opt_idx = 0;
             
             for j=1:N
                idx_p1   = find(LD(:,2)==j);    % 当前列 所有非零元素 对应的LD矩阵中位置（向量）
                idx_m = idx_p1(1);              % 当前列的第一个元素对应的序号即为D元素在LD矩阵中位置             
                i = j;
                %b_i = b_i / x_jj
                b_t(i) = b_t(i) / LD(idx_m,3);
                
                opt_idx = opt_idx + 1;
                
                idx0 = i;               % b向量列被更新元素行序号
                idx1 = 1;               % 没有使用
                idx2 = idx_m;           % b向量列更新需要使用LD矩阵中元素中序号（1维）
                idx3 = 1;               % 没有使用
                
                opt(opt_idx,1)   = 1;
                opt(opt_idx,2)   = CMD_SOLVE_DX;
                opt(opt_idx,3:6) = [idx0,idx1,idx2,idx3];
                    
                
             end
             
             x = b_t;
            
        end
        
        function [x,opt] = slove_ltx_prior(obj,N,LD,b)
            
             CMD_SOLVE_LX = 3;
             CMD_SOLVE_LTX= 4;
             CMD_SOLVE_DX = 5;

             b_t = b;               %% 待处理的右边项
             
             prior = ones(N,1);
             % 统计总的运算量
             opt_num = length(LD) - N;
             opt     = zeros(opt_num,6);
             opt_idx = 0;
             
             for i = N:-1:1         %% L转置后上三角矩阵，需要从下向上求解处理，即需要从LD矩阵的第i行向上求解
                 
                nnz_coli = length(find(LD(:,1)==i));   % 当前列非零元素数量
                idx_p1   = find(LD(:,1)==i);           % 当前列 所有非零元素 对应的LD矩阵中位置（向量）
                
                for j = 1:nnz_coli-1
                    r_idx_m =  LD(idx_p1(j),2);        % 待修正的行序号,需要考虑转置，应该是LD中的列序号
                    idx_m   = idx_p1(j);               % 待修正元素在稀疏矩阵LD中的位置            
                    
                    %b_i = b_i - b_j*x_ij
                    b_t(r_idx_m) = b_t(r_idx_m) - b_t(i)*LD(idx_m,3);
                    
                    idx0 = r_idx_m;         % b向量列被更新元素行序号
                    idx1 = i;               % b向量列更新需要使用b向量中元素的行序号
                    idx2 = idx_m;           % b向量列更新需要使用LD矩阵中元素中序号（1维）
                    idx3 = 1;               % 没有使用
                   
                    prior(r_idx_m)   = max(prior(r_idx_m),prior(i))+1;
                    
                    opt_idx          = opt_idx+1;
                    opt(opt_idx,1)   = prior(r_idx_m);
                    opt(opt_idx,2)   = CMD_SOLVE_LTX;
                    opt(opt_idx,3:6) = [idx0,idx1,idx2,idx3];  
                end
                
             end
             
             x = b_t;
            
        end
        
        
    end
end