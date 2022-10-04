%function [KUAN,L_k,bwL, pdxy]=incir_data(bws,skL)
function [cx,cy,r] = inscribed_circle_skel(bws,display)
%����bwsΪԭʼͼ��Ķ�ֵ��ͼ

tic;
%��ȡ�Ǽ�
skL = in_sketelon(bws);
%�Ǽܵ㵽��Ե���е�������������
BW = bwperim(bws); % ��ȡ�߽� ,'Canny'
%         BW=bwmorph(BW,'spur',1);
        [eg1,eg2]=find(BW==1);
        [s1,s2]=find(skL==1);
        %����Ǽ����Ե�ľ���
        ss1=(repmat(s1,[1,length(eg1)]))';
        e1=repmat(eg1,[1,length(s1)]);
        ss2=(repmat(s2,[1,length(eg1)]))';
        e2=repmat(eg2,[1,length(s1)]);
        %��С���루�뾶��
        dx=(ss1-e1).^2+(ss2-e2).^2;
       
        Mx = min(dx,[],1);%Ҷ��뾶�ı仯�������ݼ���,�����б���Ҷ����������***
        Mxs=sqrt(Mx);
     [M,w]=max(Mxs);
% if M(1)==0
%     KUAN=1;
%     %             Me=1;
% else
%     KUAN=2*max(M(1)); % *********����� Ϊʵ�������һ��
%     %             Me=Mex;
% end
% L_k=1+length(s1);

%% ��һ�ַ�ʽ
%         if KUAN>2
%             tic
% ds=samples_ypjg_dis1(bws,BW,1); % �����и��ߴ��뺯��
% %             toc
% if isempty(ds)==0
%
%     %             histogram(ds1);
%     %             xlabel('�и��߳���','FontSize',20,'FontWeight','bold');
%     %             ylabel('Ƶ��','FontSize',20,'FontWeight','bold');
%     %             legend({'ֱ����20'},'FontSize',15,'FontWeight','bold');
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

%% �ڶ��ַ�ʽ

%         if KUAN>6.8
%             ds=samples_ypjg_dis1(bws,BW,1); % �����и��ߴ��뺯��
% %             histogram(ds1);
% %             xlabel('�и��߳���','FontSize',20,'FontWeight','bold');
% %             ylabel('Ƶ��','FontSize',20,'FontWeight','bold');
% %             legend({'ֱ����20'},'FontSize',15,'FontWeight','bold');
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
%             KUAN=2*max(M(1)); % *********����� Ϊʵ�������һ��
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
        %% �Ǽ���ȡ����
        %% input:
        % bws����������Ƭ�̶�ֵͼ�񣨾����и��
        %% output:
        % skL��������ȡ�ĹǼ�
        %% ����ʼ
        sk = bwmorph(BW,'skel',Inf);  % ��ȡ�Ǽ�
        bm = ones(3); % ����3*3��ȫΪ1�ľ���
        bm(2,2) = 0;  % ���м����ֵ����Ϊ0
        [ro,la] = size(BW); % ��ֵͼ��Ĵ�С
        ze1 = zeros(ro,1);  % ����ͬ��ֵͼ��ͬ������ȫ0������
        ze2 = zeros(1,la+2); % ����һ�ж�ֵͼ����+2��������
        sk1 = [ze1,sk,ze1];  % ���Ǽ�ͼ���������߸���һ��0����
        sk2 = [ze2;sk1;ze2];  % ���Ǽ����¸���һ��0����
        sumsk = zeros(ro,la);  % ����ȫ0��ͬ��ֵͼ���С��0����
        for i = 2:ro
            for j = 2:la
                SK = sk2(i-1:i+1,j-1:j+1);  % ÿ����ȡ3*3�ľ���
                SK1 = SK .* bm;  % ���þ���������ֵֵ0
                su = sum(sum(SK1));  % ��þ���ĺ�
                if (su > 2) && (sk2(i,j) ~= 0)  % ����ʹ���2���Ҹõ�ֵ������0
                    sumsk(i-1,j-1) = su;  % �ڵ��ֵ��ͼ
                end
            end
        end
        [front1,fn1] = bwlabel(sumsk); % ͼ����
        skk = sk;  % �Ǽ�ͼ��
        szo = zeros(ro,la);  % ȫ���ֵͼ��
        for m = 1:fn1
            sk0 = skk;
            sk00 = skk;
            sk0(find(front1==m)) = 0;  % ���Ǽ�ͼ���иýڵ�ɾ��
            [front2,fn2] = bwlabel(sk0); % ���ɾ���ýڵ���ͼ��
            if fn2 < 3  % �����֧С��3������������
                continue
            end
            farea = front2(:);  % ������ת������
            [fcounts,~] = histc(farea,1:fn2);  % ͳ�Ƹ���֧�ĳ���
            [~,tn] = sort(fcounts); % �Գ��Ƚ���������С����
            %[~,st]=min(fcounts);
            for i = 1:length(tn)-2
                sk00(find(front2==tn(i))) = 0;  % �������������֧
            end
            ss = skk - sk00;  % ��֧ͼ��
            szo = ss + szo;  % ÿ���ڵ��֧
        end
        szo(find(szo ~= 0)) = 1;
        res = skk - szo;  % ȥ����֧��Ľڵ�ͼ��
        res = logical(res);
        res = bwmorph(res,'spur',1);  % ���޳�С�ķ�֧
        
    end

end
