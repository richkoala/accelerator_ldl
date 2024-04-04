clc
clear all

utils_inst      = utils();
mat_factor_inst = matrix_factor();
instr_operate_inst= instr_operate();
%% AΪN�������
N = 8;
R = sprandsym(N,0.3,0.3,1);
A = full(R);

% N = 727;
N = 557;
data_parallel = 8;      % 数据1次读取数量64bit
cmd_parallel  = 4;      % 指令1次读取数量(8+16*4)*4bit
% 


% % ԭʼPKPt������
[N,sign,eps,delta,As_rcv,ld_th]=utils_inst.load_vc_dat(N);
A    = full(sparse(As_rcv(:,1),As_rcv(:,2),As_rcv(:,3)));
% vc�����Խ����ϵ�1�������Ǵ洢ʵ�ʲ�δ����ò�����ݣ�L������������϶Խ����ϵ�Ԫ����������
% ��������Ϊld�����ʼ������


% [L,D]= ldl(A);                              % matlab ��չ LDL �ֽ����
% [Lr,Lc] = find(L~=0);A_Lv = nonzeros(L);    % ���뷽ʽ[r c v]
% L(logical(eye(size(L))))= diag(D);
% ld_th = L;

%��ȡL��������۷����������Ϣ
% nnz_idx_th=[Lr,Lc];
% A_L = tril(A);
% [As_r,As_c]= find(A_L~=0); As_v = nonzeros(A_L);
% As_rcv=[As_r,As_c,As_v];        %A����Sparse���뷽ʽ[r c v]


%��ŷֽ⺯�� LDL_symbolic()
%����:1)����ά��;2)��ֽ����rcv��ʽ3)���۷ֽ���
% nnz_idx=mat_factor_inst.ldl_symbolic(N,As_rcv,nnz_idx_th);

% utils_inst.data_dump('ld_nnz_idx_v0.txt',nnz_idx);

%LDL �ֽ���� �����㷨����ģʽ�������з����� չ������õ���LDL
% nnz_idx = ld_th(:,1:2); 
% ld = mat_factor_inst.ldl_opt(N, As_rcv, nnz_idx, ld_th);

%LDL �ֽ���� ���� ���ȼ�/������/4����ַ�ռ�
nnz_idx = ld_th(:,1:2); 
opt = mat_factor_inst.ldl_opt_prior(N, As_rcv, nnz_idx, ld_th);

%% 操作指令按照优先级进行排序
opt_prior = instr_operate_inst.order_prior(opt);

%% 指令在优先级排序的基础上按照 输出地址从大到小进行排序
opt_prior = instr_operate_inst.order_prior_addr(opt_prior,3);    %

%% 按照1bit处理 更新优先级，指令遇到优先级等于1时，需要等待所有指令写入完毕，才可以继续后续流程
opt_prior_bit1 = instr_operate_inst.instr_prior_bit1(opt_prior);

%% 填充指令
opt_align = instr_operate_inst.align(opt_prior_bit1, 4);

%% 读地址处理
data_parallel = 32;

% addr_rd_optimal(obj,opt,instr_win_h,instr_win_w, data_rd_num, adr_idx_s,adr_idx_e)
opt_adr_rd = instr_operate_inst.addr_rd_optimal(opt_align, 8, cmd_parallel, data_parallel, 3, 6);

rd_parallel = transpose(reshape(opt_adr_rd(:,7),[cmd_parallel,length(opt_adr_rd)/cmd_parallel]));

dat_rd_cnt = sum(opt_adr_rd(:,7));
disp(strcat( '并行指令个数与并行数据读取个数之比:', num2str(length(opt_adr_rd) / cmd_parallel / dat_rd_cnt)));

% addr_wr_optimal(obj,opt,instr_win_w,instr_win_h, data_wr_num, adr_idx_s)
opt_adr_wr = instr_operate_inst.addr_wr_optimal(opt_adr_rd, 8, cmd_parallel, data_parallel, 3);
dat_wr_cnt = sum(opt_adr_wr(:,9));

disp(strcat( '指令个数与数据写入次数之比:', num2str(length(opt_adr_wr) / cmd_parallel / dat_wr_cnt)));



utils_inst.data_dump('operation_code.txt',opt_prior_bit1);
utils_inst.data_dump('operation_code_fill.txt',opt_align);

