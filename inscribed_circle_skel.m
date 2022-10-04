%function [KUAN,L_k,bwL, pdxy]=incir_data(bws,skL)
function [cx,cy,r] = inscribed_circle_skel(bws,display)
%输入bws为原始图像的二值化图

tic;
%提取骨架
skL = in_sketelon(bws);
%骨架点到边缘所有点的最近距离的最大化
BW = bwperim(bws); % 提取边界 ,'Canny'
%         BW=bwmorph(BW,'spur',1);
        [eg1,eg2]=find(BW==1);
        [s1,s2]=find(skL==1);
        %计算骨架离边缘的距离
        ss1=(repmat(s1,[1,length(eg1)]))';
        e1=repmat(eg1,[1,length(s1)]);
        ss2=(repmat(s2,[1,length(eg1)]))';
        e2=repmat(eg2,[1,length(s1)]);
        %最小距离（半径）
        dx=(ss1-e1).^2+(ss2-e2).^2;
       
        Mx = min(dx,[],1);%叶面半径的变化趋势数据集合,用于判别烟叶的轮廓外形***
        Mxs=sqrt(Mx);
     [M,w]=max(Mxs);
% if M(1)==0
%     KUAN=1;
%     %             Me=1;
% else
%     KUAN=2*max(M(1)); % *********宽变量 为实际最大宽的一半
%     %             Me=Mex;
% end
% L_k=1+length(s1);

%% 第一种方式
%         if KUAN>2
%             tic
% ds=samples_ypjg_dis1(bws,BW,1); % 调用切割线代码函数
% %             toc
% if isempty(ds)==0
%
%     %             histogram(ds1);
%     %             xlabel('切割线长度','FontSize',20,'FontWeight','bold');
%     %             ylabel('频数','FontSize',20,'FontWeight','bold');
%     %             legend({'直径：20'},'FontSize',15,'FontWeight','bold');
%     %
%
%
%
%     dm=mean(ds);
%     d_s=std(ds)/dm;
%     pdxy=[dm,d_s];
% else
%     pdxy=[mean(Mxa),std(Mxa)/mean(Mxa)];
%     %                 pdxy=[KUAN,std([KUAN,L_k])/mean([KUAN,L_k])];
% end
% %         else
% %             M_dx=mean(Mxa);
% %             S_dx=std(Mxa);
% %
% %               pdxy=[M_dx,S_dx];
% %         end

%% 第二种方式

%         if KUAN>6.8
%             ds=samples_ypjg_dis1(bws,BW,1); % 调用切割线代码函数
% %             histogram(ds1);
% %             xlabel('切割线长度','FontSize',20,'FontWeight','bold');
% %             ylabel('频数','FontSize',20,'FontWeight','bold');
% %             legend({'直径：20'},'FontSize',15,'FontWeight','bold');
% %
%             dsd=sort(ds,'descend');
%             dsdsum=cumsum(dsd)/sum(dsd);
%             dsd0=sort(ds);
%             x=dsd0(2:end);
%             y=sort(dsdsum,'descend');
%             y=y(2:end);
%             Y1=log(-log(y));
%             X1=log(x);
%             X=[ones(size(Y1)) X1];
%             [b,bint,r,rint,stats] =regress(Y1,X);
%             Fds=@(x) exp((log(-log(x))-b(1))/b(2));
%
%             y05=max(y);
%             d_s= Fds(y05(1)/2);
%             xfit=0:max(x);
%             yfit=exp(-1*exp(b(1)+b(2)*log(xfit)));
%
% %             scatter(x,y,'.');
% %
% %             hold on
%             plot(xfit,yfit)
%             hold on
%             ys=0.5;
%             plot(xfit,repmat(ys,[length(xfit),1]),'--')
% %             hold on
% %
%            csdxy=d_s;
%         else
%
%             csdxy=KUAN;
%
%         end
%%
%         Mex=mean(Mxs)*2;
%         [M,w]=max(Mxs);
%         if M(1)==0
%             KUAN=1;
% %             Me=1;
%         else
%             KUAN=2*max(M(1)); % *********宽变量 为实际最大宽的一半
% %             Me=Mex;
%         end
%         L_k=1+length(s1);

cx = s2(w);
cy = s1(w);
r = M(1);
toc

% display result
if (~exist('display','var'))
    display=1;
end
if (display)
    figure
    imshow(bws)
    hold on
    nskL = bws + skL;
    nskL(nskL(:,:)==2) = 0;
    imshow(nskL)
    hold on
    plot(cx,cy,'g-*')
    hold on
    
    t = linspace(0,2*pi,60);
    xd = cx+r*cos(t);
    yd = cy+r*sin(t);
    plot(xd,yd,'r')
    hold on
    
end

    function res = in_sketelon(BW)
        %% 骨架提取程序
        %% input:
        % bws：单个区域片烟二值图像（矩形切割后）
        %% output:
        % skL：最终提取的骨架
        %% 程序开始
        sk = bwmorph(BW,'skel',Inf);  % 提取骨架
        bm = ones(3); % 生成3*3的全为1的矩阵
        bm(2,2) = 0;  % 将中间的数值设置为0
        [ro,la] = size(BW); % 二值图像的大小
        ze1 = zeros(ro,1);  % 生成同二值图相同行数的全0列向量
        ze2 = zeros(1,la+2); % 生成一行二值图列数+2的行向量
        sk1 = [ze1,sk,ze1];  % 给骨架图像左右两边各加一列0向量
        sk2 = [ze2;sk1;ze2];  % 给骨架上下各加一行0向量
        sumsk = zeros(ro,la);  % 生成全0的同二值图像大小的0矩阵
        for i = 2:ro
            for j = 2:la
                SK = sk2(i-1:i+1,j-1:j+1);  % 每次提取3*3的矩阵
                SK1 = SK .* bm;  % 将该矩阵中心数值值0
                su = sum(sum(SK1));  % 求该矩阵的和
                if (su > 2) && (sk2(i,j) ~= 0)  % 如果和大于2并且该点值不等于0
                    sumsk(i-1,j-1) = su;  % 节点二值化图
                end
            end
        end
        [front1,fn1] = bwlabel(sumsk); % 图像标记
        skk = sk;  % 骨架图像
        szo = zeros(ro,la);  % 全零二值图像
        for m = 1:fn1
            sk0 = skk;
            sk00 = skk;
            sk0(find(front1==m)) = 0;  % 将骨架图像中该节点删除
            [front2,fn2] = bwlabel(sk0); % 标记删除该节点后的图像
            if fn2 < 3  % 如果分支小于3个，程序跳出
                continue
            end
            farea = front2(:);  % 将矩阵转成向量
            [fcounts,~] = histc(farea,1:fn2);  % 统计各分支的长度
            [~,tn] = sort(fcounts); % 对长度进行排序（由小到大）
            %[~,st]=min(fcounts);
            for i = 1:length(tn)-2
                sk00(find(front2==tn(i))) = 0;  % 留下最长的两条分支
            end
            ss = skk - sk00;  % 分支图像
            szo = ss + szo;  % 每个节点分支
        end
        szo(find(szo ~= 0)) = 1;
        res = skk - szo;  % 去除分支后的节点图像
        res = logical(res);
        res = bwmorph(res,'spur',1);  % 再剔除小的分支
        
    end

end
