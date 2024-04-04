classdef instr_operate
    methods  
        %% 按照优先级进行排序处理
        function opt_pri = order_prior(obj,opt)
            [dat, idx]=sort(opt(:,1));      
            opt_pri = opt(idx,:);   
        end
        
        %% 优先级内部按照输出地址由小到大进行排序处理
        function opt_pri = order_prior_addr(obj,opt_pri,idx_addr)
            min_pri = min(opt_pri(:,1));
            max_pri = max(opt_pri(:,1));
            
            for j= min_pri:max_pri
                r_idx_pl=((find(opt_pri(:,1)==j)));               %
                [dat,r_idx_p2]= sort(opt_pri(r_idx_pl,idx_addr));        %
                t=opt_pri(r_idx_pl(r_idx_p2),:);
                opt_pri(r_idx_pl,:)= t;
            end
        end
        
        function opt_pri = instr_prior_bit1(obj,opt_pri)
            opt_num = length(opt_pri);
            for opt_idx = 1:opt_num
                if (opt_idx == 1)
                    opt_pri(opt_idx,1) =0;
                else if (opt_idx == opt_num)
                    opt_pri(opt_idx,1) =1;
                else
                    opt_pri(opt_idx,1) = opt_pri(opt_idx+1,1)-opt_pri(opt_idx,1);
                    end
                end
            end
        end

        function opt_align = align(obj, opt, align_num)
            opt_num = length(opt);
            [depth,width] = size(opt);
            na = 2^18-1;
            
            idx   = find(opt(:,1)==1);
            idx_e = [0;idx];
            
            opt_num_prior = zeros(length(idx),1);        
            for i = 1:length(idx)
                opt_num_prior(i) = idx_e(i+1) - idx_e(i);   % 每个优先级对应的指令个数
            end
            
            fill_num = mod(opt_num_prior,align_num);        % 每个优先级下的指令个数 与 align_num 求余结果 
            fill_num(find(fill_num~=0)) = align_num - fill_num(find(fill_num~=0));  % 每个优先级需要填充的指令个数
            fill_sum = sum(fill_num);                       % 需要填充的指令总数
            fill_accm = cumsum([0;fill_num]);
            idx_e_f = idx_e + fill_accm;
           
            opt_align = ones(depth+fill_sum,width)*na;   % 填充后指令，填充后均为全1值
            
            %% 完成对应的指令填充
            for i = 1:length(idx)
                opt_align(idx_e_f(i)+1 : idx_e_f(i)+opt_num_prior(i),:) = opt(idx_e(i)+1:idx_e(i)+opt_num_prior(i),:);
            end
            
        end


        function opt_addr_rd_opt = addr_rd_optimal(obj,opt,instr_win_w,instr_win_h, data_rd_num, adr_idx_s,adr_idx_e)
        
            CMD_LDL_DIV  = 1;
            CMD_MODIFY   = 2;
            CMD_SOLVE_LX = 3;
            CMD_SOLVE_LTX= 4;
            CMD_SOLVE_DX = 5;
            
            [depth,width] = size(opt);
            cmd           = opt(:,2);
            instr_win     = instr_win_w*instr_win_h;
            opt_addr_rd_opt = [opt,zeros(depth,2)];
            na            = 2^18-1;
            
            instr_win_idx = 0;
            instr_win_num = ceil(depth/instr_win);
            
            for instr_win_idx = 1 : instr_win_num
                
                if instr_win_idx == instr_win_num
                    addr_rd_win = opt((instr_win_idx-1) * instr_win+1 :                       end , adr_idx_s:adr_idx_e);
                    mode_win    = opt((instr_win_idx-1) * instr_win+1 :                       end , 2);
                else
                    addr_rd_win = opt((instr_win_idx-1) * instr_win+1 : instr_win_idx * instr_win , adr_idx_s:adr_idx_e);
                    mode_win    = opt((instr_win_idx-1) * instr_win+1 : instr_win_idx * instr_win , 2);
                end
                
                cmd_div_idx = find(mode_win==CMD_LDL_DIV);
                cmd_slx_idx = find(mode_win==CMD_SOLVE_LX);
                cmd_sdx_idx = find(mode_win==CMD_SOLVE_DX);
                cmd_stx_idx = find(mode_win==CMD_SOLVE_LTX);
                
                addr_rd_win(cmd_div_idx,3:4) = ones(length(cmd_div_idx),2).*na;
                addr_rd_win(cmd_slx_idx,  4) = ones(length(cmd_slx_idx),1).*na;
                addr_rd_win(cmd_sdx_idx,3:4) = ones(length(cmd_sdx_idx),2).*na;
                addr_rd_win(cmd_stx_idx,  4) = ones(length(cmd_stx_idx),1).*na;
                                
                addr_rd_win_h = floor(addr_rd_win / data_rd_num);
                addr_rd_win_h_unique = unique(addr_rd_win_h);
                addr_rd_win_h_unique(find(addr_rd_win_h_unique==(floor(na/ data_rd_num)) ))=[];    %% 剔除非法数据
                    
                for i = 1: length(addr_rd_win_h_unique)
                    [idx_r,idx_c]    = find( addr_rd_win_h == addr_rd_win_h_unique(i) );
                    if instr_win_idx==9 && i==12
                        instr_win_idx
                    end
                    %指令并行度为1的指令行序号
                    %指令并行度为1的指令下地址列序号
                    if (length(idx_r)==1)
                        opt_addr_idx = [idx_r,idx_c];
                    else
                        opt_addr_idx = max([idx_r,idx_c]);
                    end
                    opt_addr_idx_r = opt_addr_idx(1);
                    opt_addr_idx_c = opt_addr_idx(2);

                    idx_rd = (instr_win_idx-1) * instr_win + opt_addr_idx_r;
                    opt_addr_rd_opt(idx_rd,7) = 1;
                    opt_addr_rd_opt(idx_rd,8) = opt_addr_idx_c;
                end
            end
        end
        
        % 在opt最后一列补1，表明为需要最后可以并行写入的地址
        function opt_addr_wr_opt = addr_wr_optimal(obj,opt,instr_win_w,instr_win_h, data_wr_num, adr_idx_s)
        
            CMD_LDL_DIV  = 1;
            CMD_MODIFY   = 2;
            CMD_SOLVE_LX = 3;
            CMD_SOLVE_LTX= 4;
            CMD_SOLVE_DX = 5;
            
            [depth,width] = size(opt);
            cmd           = opt(:,2);
            instr_win     = instr_win_w*instr_win_h;
            opt_addr_wr_opt = [opt,zeros(depth,1)];
            na            = 2^18-1;
            instr_win_idx = 0;
            instr_win_num = ceil(depth/instr_win);
            
            for instr_win_idx = 1 : instr_win_num
                
                if instr_win_idx == instr_win_num
                    addr_wr_win = opt((instr_win_idx-1) * instr_win+1 :                       end , adr_idx_s);
                else
                    addr_wr_win = opt((instr_win_idx-1) * instr_win+1 : instr_win_idx * instr_win , adr_idx_s);
                end

                addr_wr_win_h = ceil(addr_wr_win / data_wr_num)-1;
                addr_wr_win_h_unique = unique(addr_wr_win_h);

                for i = 1: length(addr_wr_win_h_unique)
                    [idx,dat]    = find(addr_wr_win_h==addr_wr_win_h_unique(i));
                    addr_wr_flag = max(idx);
                    opt_addr_idx_r = addr_wr_flag;          %指令并行度为1的指令行序号

                    idx_wr = (instr_win_idx-1) * instr_win + opt_addr_idx_r;
                    opt_addr_wr_opt(idx_wr,end) = 1;
                end
                
            end
        end
        
        
    end
end