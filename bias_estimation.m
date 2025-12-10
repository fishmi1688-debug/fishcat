function preproc=bias_estimation(rhip,lhip)


% load('lsl_itf4.mat')
% %
% t=itf4.data.t;
% lhip=itf4.data.sang(:,1);
% rhip=itf4.data.sang(:,3);
% footP=itf4.data.footRT;
% Euler=itf4.data.Euler;
% lknee=itf4.data.sang(:,2);
% nonparetic side joint minus paretic side joint
q=rhip-lhip;           %关节角度差值
len=length(q);         %数据长度
% q filter
sys=tf(1,[0.1 1]);      %滤波器传递函数 
coef=c2d(sys,0.01,'tutin');        %离散化
qfil.num=coef.Numerator{1};        %滤波器分子系数
qfil.den=coef.Denominator{1};      %滤波器分母系数
% dq filter
sys=tf([1 0],[0.05 1]);             %滤波器传递函数
coef=c2d(sys,0.01,'tutin');        %离散化
dqfil.num=coef.Numerator{1};       %滤波器分子系数
dqfil.den=coef.Denominator{1};     %滤波器分母系数
% qb filter
sys=tf(1,[0.1 1]);                  %滤波器传递函数
coef=c2d(sys,0.01,'tutin');        %离散化
qbfil.num=coef.Numerator{1};       %滤波器分子系数
qbfil.den=coef.Denominator{1};     %滤波器分母系数

%% estimate q bias

% q filter
qf=zeros(len,1);          %滤波器输出初始化
for i=1:len               %滤波器迭代
    qf(i)=qEstimation(q(i),qfil);        %滤波器输出
end
% dq filter
dq=zeros(len,1);             %滤波器输出初始化
for i=1:len                  %滤波器迭代
    dq(i)=dqEstimation(qf(i),dqfil);          %dq滤波器输出
end

% 50*10*0.01=5s 
q_resample=zeros(30,1);             %采样缓存初始化
count=0;                            %采样计数初始化
bias=0;                             %偏置初始化
q_max=0;                            %最大值初始化
q_min=0;                            %最小值初始化
q_amp=0;                            %幅值初始化
qbias=zeros(len,1);                 %偏置缓存初始化
qmax=zeros(len,1);                  %最大值缓存初始化
qmin=zeros(len,1);                  %最小值缓存初始化
qamp=zeros(len,1);                  %幅值缓存初始化
t=0:0.01:(len-1)*0.01;              %时间向量
for i=1:len             %偏置估计迭代
    count=count+1;           %采样计数加1
    if(count==10)       %每10个采样点更新一次
        count=0;        %采样计数清零
        % resample
        q_resample=circshift(q_resample,1);     %缓存循环移位
        q_resample(1)=qf(i);           %最新采样值存入缓存
        q_max=max(q_resample);      %计算最大值
        q_min=min(q_resample);      %计算最小值
        bias=0.5*(q_max+q_min);     %计算偏置
        q_amp=0.5*(q_max-q_min);    %计算幅值
    end
    qbias(i) = biasEstimation(bias,qbfil);      %偏置滤波器输出
    qmax(i)=q_max;             %最大值存储
    qmin(i)=q_min;       %最小值存储
    qamp(i)=q_amp;       %幅值存储  
end

qf_nb=qf-qbias;     %去偏滤波器输出

preproc.qf=qf;      %滤波器输出
preproc.qf_nb=qf_nb; %去偏滤波器输出
preproc.qmax=qmax;   %最大值
preproc.qmin=qmin;   %最小值
preproc.dq=dq;       %dq滤波器输出
preproc.qamp=qamp;   %幅值
figure              
x1=subplot(3,1,1);                  %绘图
hold on                     
plot(t,q)
plot(t,qf)
plot(t,preproc.qamp)
legend('q','qf','amp')
x2=subplot(3,1,2);
hold on
plot(t,gradient(qf)*100)
plot(t,dq)
legend('dq','dqf')
x3=subplot(3,1,3);
plot(t,preproc.qf_nb)
hold on
plot(t,preproc.qf)
plot(t,preproc.qmax)
plot(t,preproc.qmin)
legend('q no bias', 'q bias','q max','q min')
linkaxes([x1,x2,x3],'x')
%% 



end



function y=qEstimation(u,qfil)          %滤波器实现
persistent u_last y_last                %持久变量
if isempty(u_last)              %初始化 
    u_last=0;                   
end
if isempty(y_last)               %初始化    
    y_last=0;                   %初始化
end

% y=0.9802*y_last+0.009901*u_last+0.009901*u;
y=-qfil.den(2)*y_last+qfil.num(2)*u_last+qfil.num(1)*u;       %滤波器计算  
u_last=u;       %更新持久变量
y_last=y;       %更新持久变量
end

function y=dqEstimation(u,dqfil)        %滤波器实现
persistent u_last y_last        %持久变量
if isempty(u_last)        %初始化
    u_last=0;           
end
if isempty(y_last)      %初始化 
    y_last=0;
end
y=-dqfil.den(2)*y_last+dqfil.num(2)*u_last+dqfil.num(1)*u;   %滤波器计算
u_last=u;       %更新持久变量
y_last=y;       %更新持久变量
end

function y=biasEstimation(u,qbfil)     %滤波器实现
persistent u_last y_last    %持久变量
if isempty(u_last)      %初始化 
    u_last=0;
end
if isempty(y_last)      %初始化
    y_last=0;
end
y=-qbfil.den(2)*y_last+qbfil.num(2)*u_last+qbfil.num(1)*u;   %滤波器计算
u_last=u;       %更新持久变量
y_last=y;       %更新持久变量
end



