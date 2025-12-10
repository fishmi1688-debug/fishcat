function footp=GaitPExtraction(t,footF,gflag)    %脚底压力步态相位提取

% load('lsl_itf4.mat')
% t=itf4.data.t;
% lhip=itf4.data.sang(:,1);
% rhip=itf4.data.sang(:,3);
% footF=itf4.data.footRT;
% Euler=itf4.data.Euler;
% lknee=itf4.data.sang(:,2);

%
footT=150;          % 足底压力阈值
timeT=0.25;         % 时间阈值
len=length(footF);         % 数据长度
footP_logic=zeros(len,1);   %足底压力逻辑变量
% obtain swing phase and stand phase
time_c=0;     %电流时间记录
phase_c=0; % 当前相位记录 0 stance 1 swing
pha=zeros(len,1);     %相位记录 0 stance 1 swing
w=pi;          %初始频率
w_init=2*pi;        %初始频率
fre=zeros(len,1); %     频率记录
freg=w_init/2*ones(len,1); %    估计频率记录
stri=zeros(len,1);          %步数记录
for i=1:len                 %数据遍历
    
    if(gflag(i)==1)         % 行走状态
    if(footF(i)<footT)         % 脚底压力小于阈值
        % swing
        footP_logic(i)=1;           %摆动
    else
        % stance
        footP_logic(i)=0;        %支撑  
    end
        % phase switch 
    if(i>=2)
        if(footP_logic(i)~=footP_logic(i-1))        %相位切换
            time_l=time_c;       %上次切换时间记录
            time_c=t(i);     % 本次切换时间记录
            %   步态相位更新  
            if((time_c-time_l)<timeT) 
                %   时间未达到阈值
            else
                %   时间达到阈值
                phase_c=footP_logic(i);       %相位更新           
            end
        end
    end
    else
        %   非行走状态
        phase_c=0;     %相位归零
    end
    
    pha(i)=phase_c;     %相位记录
end

%% obtain w w_g
step=0;         %步数记录
time_itv=zeros(3,1);     %步态时间间隔记录
for i=2:len
    if(gflag(i)==1)      % 行走状态
        % obtain the step information
        if(pha(i)~=pha(i-1))     %相位切换
            step=step+1;
            time_itv(1)=time_itv(2);     %时间间隔更新
            time_itv(2)=time_itv(3);     %时间间隔更新
            time_itv(3)= t(i);     %时间间隔更新
            
            %   频率更新
            if(step<=2)         %步数小于等于2
            w=2*pi/(time_itv(3)-time_itv(2));       %频率更新
            w_g=w_init/2;       %估计频率更新
            else
            w=2*pi/(time_itv(3)-time_itv(2));       %频率更新
            w_g=2*pi/(time_itv(3)-time_itv(1));     %估计频率更新
            end
        
        end

        
    else
        step=0;     %步数归零
        w=w_init;    %频率归初始值
        w_g=w_init/2;       %估计频率归初始值
    end
    
    
    fre(i)=w;        %频率记录
    freg(i)=w_g;     %估计频率记录
    stri(i)=step;    %步数记录
end


%% obtain phase p   
dt=0.01;            %采样时间间隔
p=zeros(len,1);     %相位记录
for i=2:len               %数据遍历
    if(gflag(i)==0)         %非行走状态
        p(i)=0;             %相位归零
        w_g=0;              %频率归零   
    else
        % switch from stance to swing
        w_g=freg(i);            %频率更新
    end
    
    if(pha(i-1)==0&&pha(i)==1)      %   从支撑到摆动
    p(i)=0;             %相位归零
    else
        p(i)=p(i-1)+w_g*dt;     %相位更新
        if(p(i)>2*pi)           %相位归一化
           p(i)=2*pi;           %相位归一化
        end
    end
end

%% phase correction
xi=3000;            %相位校正参数   
coef=1./(1+(footF-footT).^2/xi);            %相位校正系数


footp.stri=stri;        %步数
footp.fre=fre;      %频率
footp.freg=freg;        %估计频率
footp.pha=pha;          %相位
footp.p=p;              %相位
footp.footT=footT;      %足底压力阈值
footp.coef=coef;        %相位校正系数
%%


p_c=zeros(len,1);           %校正后相位记录

for i=2:len          %数据遍历
    if(gflag(i)==0)           %非行走状态
     p_c(i)=footp.p(i);         %相位不变
    else
        if(footp.p(i)>0.75*2*pi&&footp.coef(i)>0.05)        %相位大于0.75*2*pi且校正系数大于0.05
            lambda=(footp.coef(i));         %校正系数
            p_c(i)=lambda*(2*pi-0.314*(1-lambda))+(1-lambda)*footp.p(i);        %相位校正
            if(p_c(i)<p_c(i-1))      %相位单调性保持
            p_c(i)=p_c(i-1);        %相位单调性保持
            end
        else
            p_c(i)=footp.p(i);      %相位不变
        end
    end
    
end


footp.p_c=p_c;          %校正后相位



figure
x1=subplot(2,1,1);
hold on
plot(footF/600)
plot(footT*ones(len,1)/600)
plot(footP_logic)
legend('foot p','The','phase naive')
x2=subplot(2,1,2);
hold on
plot(pha)
plot(freg/2/pi)
plot(footp.p/2/pi)
legend('phase','wg','p')
linkaxes([x1,x2],'x')



figure
hold on
plot(t,footp.p/2/pi,'black','Linewidth',0.5)
plot(t,footp.p_c/2/pi,'red','LineWidth',1)
plot(t,footp.coef,'blue','LineWidth',0.3)
legend('$foot_p/(2\pi)$','$foot_{pc}/(2\pi)$','$\lambda$','interpreter','latex')
xlabel('time (s)', 'interpreter','latex')
set(gca,'fontSize',16)

end