function [r,cx,cy] = inscribed_circle_opti(contourImage,display)
% --- ����˵�� ---
% ��ȡ����ͼ�ε��������Բ
% �������Թ滮���⣬Ŀ�꺯��Ϊ��max r(���뾶)
% Լ��������min(sqrt((x-a)^2+(y-b)^2))>=r
% ����(x,y)Ϊͼ�α߽��ϵ�����꣬(a,b)ΪԲ�����꣬rΪԲ�뾶��
%
% --- ���� ---
% ContourImage : Ŀ��ͼ������ͼ
% display : �Ƿ����ͼ�ν��չʾ
%
% --- ��� ---
% R : �������Բ�뾶 
% cx : ����ԲԲ��x���� 
% cy : ����ԲԲ��y���� 
%
% Programmer: Robin An, 2021-10-01
% last modified by Robin An on 2021-10-01
% reference: http://www.matlabsky.com/thread-39620-1-1.html
%            ���Ƽ���ѧ��ѧϵ л�л�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fun = @(x)-x(3); %x = beta = [beta(1),beta(2),beta(3)]
%BW = imread(filename);
[y,x] = find(contourImage);
warning off

% fun Ŀ�꺯��
% funcontestԼ������
% fmincon(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
initVal = [mean(x),mean(y),50];
lb = [min(x);min(y);0];
ub = [max(x);max(y);min(max(y)-min(y),max(x)-min(x))];%���������ֵ������Сһ��
[beta,fval] = fmincon(fun,initVal,[],[],[],[],lb,ub,@funcontest);
% fun������Ϊbeta��beta(1)��beta(2)��beta(3)��

cx = beta(1);
cy = beta(2);
r = beta(3);
toc

% res = inpoly([cx,cy],[x,y])
% res1 = inpoly_se(cx,cy,x,y)
% res2 = inpoly_bl(cx,cy,x,y)
% res3 = inpoly([406,384],[x,y])
% res4 = inpoly_se(406,384,x,y)
% res5 = inpoly_bl(406,384,x,y)
% figure
% imshow(contourImage)
% hold on
% plot(cx,cy,'g-*') % 282,236
% hold on
% plot(406,384,'b-+')

% display result
if (~exist('display','var'))
    display=1;
end
if (display)
    figure
    imshow(contourImage)
    hold on
    plot(beta(1),beta(2),'g-*')
    hold on
    
    t = linspace(0,2*pi,60);
    xd = beta(1)+beta(3)*cos(t);
    yd = beta(2)+beta(3)*sin(t);
    plot(xd,yd,'r')
    hold off
end

% ������Լ���Ӻ�����C���
    function [C,Ceq] = funcontest(beta)
        % a = beta(1),b = beta(2),r = beta(3)
                
        C = beta(3)-min(sqrt((beta(1)-x).^2+(beta(2)-y).^2));
        
        % ���������ڲ�
        %Ceq = [];
        Ceq = inpoly([beta(1),beta(2)],[x,y]) - 1;
    end
end
