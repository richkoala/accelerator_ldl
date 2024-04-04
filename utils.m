classdef utils
    methods
        %% 1�� ����solve���Է����������������L/D������� lx/dx/ltx
        function [N,sign,eps,delta,As,LDs]=load_vc_dat(obj,dim)
            
            PKPt   = load(strcat('vc2mat_',num2str(dim),'/PKPt.txt'));
            sign   = load(strcat('vc2mat_',num2str(dim),'/Sign.txt'));
            eps    = load(strcat('vc2mat_',num2str(dim),'/eps.txt'));
            delta  = load(strcat('vc2mat_',num2str(dim),'/delta.txt'));
            L      = load(strcat('vc2mat_',num2str(dim),'/L.txt'));
            D      = load(strcat('vc2mat_',num2str(dim),'/D.txt'));
            
            N    = length(sign);
            As_c = PKPt(1:end-1,1);   %����PKPtΪ�����Ǿ���ķ���Ԫ�أ���Ҫ����ת�ô���col����Ϊrow
            As_r = PKPt(1:end-1,2);   %����PKPtΪ�����Ǿ���ķ���Ԫ�أ���Ҫ����ת�ô���row����Ϊcol
            As_v = PKPt(1:end-1,3);
            As   = [As_r,As_c,As_v];
 
            Lr = L(1:end-1,1);
            Lc = L(1:end-1,2);
            Lv = L(1:end-1,3);
            
            LDr = [[1:N]';Lr];
            LDc = [[1:N]';Lc];
            LDv = [D     ;Lv];    
            [dat,idx] = sort(LDc); % ������������D�����ݲ���
            LDr = LDr(idx);
            LDc = LDc(idx);
            LDv = LDv(idx);
            LDs = [LDr,LDc,LDv];
 
        end
           
        %% 2�� ��������
        function data_dump(obj,file_name,dat)
            [depth,width] = size(dat);
            fid = fopen(file_name,'wt');
            for i = 1:depth
                for j = 1:width
                    fprintf(fid,'%d\t',dat(i,j));
                    if (j==width)
                        fprintf(fid,'\n');
                    end
                end
            end
            fclose all;
        end
        
    end
end
            