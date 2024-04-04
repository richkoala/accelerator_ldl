classdef matrix_factor
    methods  
        %% 符号分解函数 LDL_symbolic(),符号分解函数可能与LDL分解的存在差异，需要注意
        function nnz_idx=ldl_symbolic(obj,N,As_rcv,nnz_idx_th)
            %获取A下三角矩阵A L的稀疏编码
            nnz_A_L = length(As_rcv);   
            nnz_idx  = As_rcv(:,1:2);    %所有非零点的行列坐标
            nnz_idx_add =[];             %新增计算点位
            for j= 1:N
                if(mod(j,2)==0)
                    disp(sprintf('符号分解 当前开展第 %d 列',j));
                end
                [dat,idx] = sort(nnz_idx(:, 2)); %按照列方向重新从小到大获取序号
                nnz_idx_t = nnz_idx(idx,:);      %列方向重新排序处理
                nnz_idx   = nnz_idx_t;
                
                r_idx_pl = (find(nnz_idx(:,2)==j));     %   获取当前列的所有行序号
                [dat,r_idx_p2] = sort(nnz_idx(r_idx_pl,1));  %按照行序号进行从小到大排序处理，获取行序号重排所对应的序号
                t = nnz_idx(r_idx_pl(r_idx_p2),:);    
                nnz_idx(r_idx_pl,:) = t;

                nnz_colj = length(find(nnz_idx(:,2)==j));    %获取当前列非零元素长度
                r_colj   = nnz_idx((find(nnz_idx(:,2)==j)),1);  %获取当前列的行序号

                for i = 2 : nnz_colj %% 第一个元素为对角线元素，可以忽略
                    for k = i+1 : nnz_colj
                        
                        A_L_midx=[r_colj(k),r_colj(i)];
                        if(find(ismember(nnz_idx,A_L_midx,'rows'))>0) %%没有对应元素时进行添加处理
                        else
                            nnz_idx    = [nnz_idx;    A_L_midx];
                            nnz_idx_add= [nnz_idx_add;A_L_midx];
                        end
                    end
                end
            end
            % 测试验证正确性
            disp(sprintf('step1:符号分解处理，判断新增元素的正确性,偏差=%d' ,max(max(abs (nnz_idx_th-nnz_idx)))));
        end
        
        %% LDL 分解计算右视算法处理模式，按照列方向处理 展开计算得到的LDL
        function ld= ldl_opt (obj,N,As_rcv,nnz_idx, ld_th)
        %1) 初始化时值损作。
            nnz_ld =length(nnz_idx);
            nnz_A_L=length(As_rcv);
            A_Ev=zeros(nnz_ld,1);   %% 扩展后的A矩阵元素
            for i = 1:nnz_A_L       %% 外层--循环1)待分解的列
                idx_ini = find(ismember(nnz_idx, [As_rcv(i,1),As_rcv(i,2)],'rows'));
                A_Ev(idx_ini)= As_rcv(i,3);
            end
            A_Es = [nnz_idx,A_Ev];
        %2.1)稀疏矩阵的LDL计算
            for j= 1:N      %%外层--循环1)待分解的列
                if(mod(j,20)==0)
                    disp(sprintf('当前开展第 %d 列更新',j));
                end
                nnz_colj = length(find(A_Es(:,2)==j));   %当前列非零元素数量
                idx_p1   = find(A_Es(:,2)==j);           %当前列非零元素对应的行序号

                for i=2:nnz_colj    %%行方向 (行找到非零元素)
                    idx_m =idx_p1(i);
                    idx_l0 =idx_p1(1);
                    A_Es(idx_m,3) = A_Es(idx_m,3)/A_Es(idx_l0,3);
                end

                for jk= 2:nnz_colj
                    for ik=jk:nnz_colj
                        r_idx_m = A_Es(idx_p1(ik),1);
                        c_idx_m = A_Es(idx_p1(jk),1);
                        idx_l0  = idx_p1(1);
                        idx_l1  = find(ismember(nnz_idx,[c_idx_m,j],'rows'));
                        idx_12  = find(ismember(nnz_idx,[r_idx_m,j],'rows'));
                        idx_m   = find(ismember(nnz_idx,[r_idx_m,c_idx_m],'rows'));

                        A_Es(idx_m,3) = A_Es(idx_m,3) - A_Es(idx_l1,3) * A_Es(idx_12,3) * A_Es(idx_l0,3);
                    end
                end
            end
            ld  = A_Es;
%             ld_dense = full(sparse(ld(:,1),ld(:,2),ld(:,3)));
            disp(sprintf('step2.1;[DL分解处理，判断分解后LD元素计算偏差=%e',max(max(abs(ld_th(:,3)-ld(:,3))))));
        end
        
        %% LDL分解引入优先级
        function opt = ldl_opt_prior(obj,N,As_rcv,nnz_idx, ld_th)
        %1) 初始化时值损作。
            nnz_ld =length(nnz_idx);
            nnz_A_L=length(As_rcv);
            A_Ev=zeros(nnz_ld,1);   %% 扩展后的A矩阵元素
            for i = 1:nnz_A_L       %% 外层--循环1)待分解的列
                idx_ini = find(ismember(nnz_idx, [As_rcv(i,1),As_rcv(i,2)],'rows'));
                A_Ev(idx_ini)= As_rcv(i,3);
            end
            A_Es = [nnz_idx,A_Ev];
            CMD_LDL_DIV = 1;
            CMD_MODIFY  = 2;
            opt_num     = 0;
            
            %统计总的运算量
            for j=1:N
                nnz_colj = length(find(A_Es(:,2)==j));
                opt_num  = opt_num + nnz_colj - 1;
                opt_num  = opt_num + nnz_colj*(nnz_colj-1)/2;
            end
            
            prior   = ones(length(A_Es),1);
            opt     = zeros(opt_num,6);
            opt_idx = 0;
            
        %2.1)稀疏矩阵的LDL计算
            for j= 1:N      %%外层--循环1)待分解的列
                if(mod(j,20)==0)
                    disp(sprintf('当前开展第 %d 列更新',j));
                end
                
                nnz_colj = length(find(A_Es(:,2)==j));   % 当前列非零元素数量
                idx_p1   = find(A_Es(:,2)==j);           % 当前列非零元素对应的行序号

                for i=2:nnz_colj                         % 行方向 (行找到非零元素)
                    opt_idx = opt_idx + 1;  
                    idx_m   = idx_p1(i);
                    idx_l0  = idx_p1(1);
                    A_Es(idx_m,3) = A_Es(idx_m,3)/A_Es(idx_l0,3);
                    
                    prior(idx_m)   = max(prior(idx_m),prior(idx_l0))+1; % 优先级设计
                    opt(opt_idx,1) = prior(idx_m);
                    opt(opt_idx,2) = CMD_LDL_DIV;
%                     opt(opt_idx,3:6)=[idx_m,idx_l0,sign_dat,65535,idx_m];
                    opt(opt_idx,3:6)=[idx_m,idx_l0,0,0];    %地址编码
                    
                end

                for jk= 2:nnz_colj
                    for ik=jk:nnz_colj
                        opt_idx = opt_idx + 1;
                        r_idx_m = A_Es(idx_p1(ik),1);
                        c_idx_m = A_Es(idx_p1(jk),1);
                        idx_l0  = idx_p1(1);
                        idx_l1  = find(ismember(nnz_idx,[c_idx_m,j],'rows'));
                        idx_12  = find(ismember(nnz_idx,[r_idx_m,j],'rows'));
                        idx_m   = find(ismember(nnz_idx,[r_idx_m,c_idx_m],'rows'));

                        A_Es(idx_m,3) = A_Es(idx_m,3) - A_Es(idx_l1,3) * A_Es(idx_12,3) * A_Es(idx_l0,3);
                        
                        prior(idx_m)   = max([prior(idx_m),prior(idx_l1),prior(idx_12),prior(idx_l0)])+1; %优先级设计
                        opt(opt_idx,1) = prior(idx_m);
                        opt(opt_idx,2) = CMD_MODIFY;
                        opt(opt_idx,3:6) =[idx_m,idx_l0,idx_l1,idx_12]; %地址表为[待修正元素地址/D元素/L0元素/L1元素]
                        
                    end
                end
            end
            ld  = A_Es;
%             ld_dense = full(sparse(ld(:,1),ld(:,2),ld(:,3)));
            disp(sprintf('step2.2;[DL分解处理，判断分解后LD元素计算偏差=%e',max(max(abs(ld_th-ld)))));
        end
        
    end
end
            