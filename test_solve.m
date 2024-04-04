clc
clear all

mat_solver_inst = matrix_slover();
utils_inst = utils();
%% A为N阶正定矩阵
N = 100;

R = sprandsym(N,0.3,0.3,1);
A = full(R);
x = rand(N,1);


[L,D] = ldl(A);
b_lx   = L*x;
b_ldx  = L*D*x;
b_ldlx = L*D*L'*x;

[Lr,Lc] = find(L~=0);
Lv = nonzeros(L);

[Lr,Lc] = find((L+D)~=0);
L(logical(eye(size(L))))= diag(D);  % 将对角线1替换为D元素
Lv      = nonzeros(L);
LD      = [Lr,Lc,Lv];

[x_l,opt] = mat_solver_inst.slove_lx_prior(N,LD,b_lx);
max(abs(x-x_l))

[x_l,opt_lx] = mat_solver_inst.slove_lx_prior(N,LD,b_ldx);
[x_d,opt_dx] = mat_solver_inst.slove_dx_prior(N,LD,x_l);
max(abs(x_d-x))

[x_l  , opt_lx] = mat_solver_inst.slove_lx_prior(N,LD,b_ldlx);
[x_d  , opt_dx] = mat_solver_inst.slove_dx_prior(N,LD,x_l);
[x_lt ,opt_ltx] = mat_solver_inst.slove_ltx_prior(N,LD,x_d);
max(abs(x_lt-x))


%% 按照指令优先级进行排序
[dat, idx]=sort(opt_lx(:,1));      
opt_lx_1 = opt_lx(idx,:);             

%% 指令优先级内部按照被更新元素位置进行排序 
min_pri = min(opt_lx_1(idx,1));
max_pri = max(opt_lx_1(idx,1));
for j= min_pri:max_pri
    r_idx_pl=((find(opt_lx_1(:,1)==j)));               %获取当前列的所有行序号
    [dat,r_idx_p2]= sort(opt_lx_1(r_idx_pl,3));        %按照行序号进行从小到大排序处理，获取行序号重排所对应的序号
    t=opt_lx_1(r_idx_pl(r_idx_p2),:);
    opt_lx_1(r_idx_pl,:)= t;
end






