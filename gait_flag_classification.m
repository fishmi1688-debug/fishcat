function gflag=gait_flag_classification(t,preproc,footF)        %步态开始停止分类


%% gait start and gait stop classifier
qf_nb=preproc.qf_nb;        %滤波后膝关节角度
dq=preproc.dq;              %滤波后膝关节角速度
%% 1: r obtainment
len=length(qf_nb);          %数据长度
r=zeros(len,1);       %r初始化
for i=1:len              %r计算
    r(i)=sqrt(qf_nb(i)^2+0.1*dq(i)^2);          %r计算公式
end

%% 2 logic
% it is classified by a fuzzy logic:
%%gait start
% start r is bigger than a certain level
% gait switch from stance to swing or gait swing time bigger than a level
%%gait stop 
% time duration of stance phase is bigger than a certain level
% r is continuing less than a certain level
footTH=150;         %脚底压力阈值
gflag=zeros(len,1);         %步态标志初始化
g_flag=0;                   %当前步态状态
rT=0.3;                     %r阈值
gaitSD=1.0;                 %步态持续时间阈值
timer_r=0;                  %r计时器
time_foot=0;                %脚底压力计时器
for i=2:len        %步态分类循环
    % start detection
    if(gflag(i-1)==0)        %当前步态为停止
       if(r(i)>rT)               %r大于阈值
          timer_r=0;             %r计时器清零
       else
          timer_r=timer_r+0.01;     %r计时器累加
       end
       % switch from stance to swing
       if(footF(i-1)<=footTH&&footF(i)<=footTH)         %脚底压力小于阈值
           time_foot=time_foot+0.01;                      %脚底压力计时器累加
       else
           time_foot=0;                                  %脚底压力计时器清零
       end
       % gait start classify
       if((footF(i-1)>footTH&&footF(i)<=footTH))            %脚底压力由大于阈值变为小于阈值
           ggflag=1;
       end
       if((footF(i-1)>footTH&&footF(i)<=footTH&&timer_r<2))      %脚底压力由大于阈值变为小于阈值且r计时器小于2s
           % walking
           g_flag=1;                   %步态状态设为行走
           time_foot=0;               %脚底压力计时器清零
           timer_r=0;                 %r计时器清零
       end       
       
    end
    
    
    % stop detection
    if (gflag(i-1)==1)              %当前步态为行走
        
       if(r(i)>rT)               %r大于阈值
          timer_r=0;             %r计时器清零
       else
          timer_r=timer_r+0.01;     %r计时器累加
       end
       
       if(footF(i-1)>footTH&&footF(i)>footTH)
           time_foot=time_foot+0.01;                      %脚底压力计时器累加
       else
           time_foot=0;                                  %脚底压力计时器清零
       end
       
       % gait stop classify
       if(timer_r>2.5||(time_foot>2.0&&timer_r>2.0))    %r计时器大于2.5s或脚底压力计时器大于2s且r计时器大于2s   
           % walking
           g_flag=0;            %步态状态设为停止
           time_foot=0;         %脚底压力计时器清零
           timer_r=0;           %r计时器清零    
       end     
       
    end
    
    gflag(i)=g_flag;        %步态状态赋值
end

figure
hold on
plot(t,footF/300)                          %脚底压力归一化绘图
plot(t,footTH*ones(len,1)/300)             %脚底压力阈值归一化绘图
plot(t,r)                                  %r绘图
plot(t,gflag)                              %步态标志绘图
legend('footF/300','T/300','r','gflag')             


end