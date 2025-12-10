function gait=adaptive_oscillators(preproc,gflag,footp,Euler)           %步态相位估计

qf_nb=preproc.qf_nb;        %去噪关节角度
len=length(qf_nb);          %数据长度

dt=0.01;                    %采样时间
t=0:dt:dt*(len-1);          %时间序列
%
order=2;                    %振荡器阶数
ao.para=zeros(order,3);     %振荡器参数初始化 [amp freq phase]
q_esti=zeros(1,order);          %估计关节角度初始化
% 
gait.p_pure=zeros(len,1);       %步态相位纯值
gait.p=zeros(len,1);            %步态相位mod 2pi
gait.w=zeros(len,1);            %步态频率
gait.q=zeros(len,1);            %估计关节角度
gait.a0=zeros(len,1);           %振荡器0阶幅值
for i=1:len                     
    %% adaptive oscillators
    %init
    if(gflag(i)==1&&gflag(i-1)==0&&i>1)         %   启动时初始化
    q= preproc.qf_nb (i);           %当前去噪关节角度
    dq=preproc.dq(i);           %当前关节角速度
    %
    T=2*pi/footp.freg(i);           %当前步态周期
    T_aefa=2*T;                 %aefa时间常数
    T_amig=2*T;                 %amig时间常数
    ao.yita=2/T_aefa;           %aefa收敛速度
    ao.vaomig=20/(T_amig)^2;    %amig收敛速度
    ao.vfai=sqrt(24.2*ao.vaomig); %vfai参数 设置
    % init
    ao.para(1,1)=0;             %   0阶振荡器幅值设为0
    ao.para(1,3)=pi/2;          %   0阶振荡器相位设为pi/2
    ao.para(2,1)=0.4;               %   1阶振荡器幅值初始设为0.4
    ao.para(2,2)=footp.freg(i);             %   1阶振荡器频率初始设为当前步态频率
    ao.para(2,3)=atan2(ao.para(2,2)*q,dq);      %   1阶振荡器相位初始设为atan2(omega*q,dq)
    % correct the amp
    ao.para(2,1) = max(q/sin(ao.para(2,3)),0.3);        %   1阶振荡器幅值设为q/sin(fai)
    end
            
    if(gflag(i)==1)
            % tracking
            q= qf_nb (i);           %当前去噪关节角度   
            for j=1:order
            q_esti(j)=ao.para(j,1)*sin(ao.para(j,3));       %估计关节角度
            end 
            ao.y=sum(q_esti);       %振荡器输出
            F=q-ao.y;                   %估计误差
            % aefa_0_dot
            ao.para_dot(1,1)=ao.yita*F;         %0阶振荡器幅值变化率
            %aerfa_dot
            for j=2:order
            ao.para_dot(j,1)=ao.yita*F*sin(ao.para(j,3));           %j阶振荡器幅值变化率
            end
            %aomig_dot
            %para_dot(2,2)=vaomig*F/sum(para(:,1))*cos(para(2,3));
            ao.para_dot(2,2)=ao.vaomig*F*cos(ao.para(2,3));         %2阶振荡器频率变化率
            % fai_dot
            for j=2:order
            %para_dot(j,3)=para(j,2)+vfai*F/sum(para(:,1))*cos(para(j,3)); 
            ao.para_dot(j,3)=ao.para(j,2)+ao.vfai*F*cos(ao.para(j,3));          %j阶振荡器相位变化率
            end
            for j=2:order
            ao.para_dot(j,3)=ao.para(j,2)+ao.vfai*F*cos(ao.para(j,3));          %j阶振荡器相位变化率
            if(ao.para_dot(j,3)<0)
                ao.para_dot(j,3)=0;              %相位变化率非负限制
            end
            end
            % update para
            %alpha_0
            ao.para(1,1)=ao.para(1,1)+ao.para_dot(1,1)*dt;              %0阶振荡器幅值更新
            ao.para(1,1)=0;                             %0阶振荡器幅值设为0
                %alpha
            for j=2:order
            ao.para(j,1)=ao.para(j,1)+ao.para_dot(j,1)*dt;          %j阶振荡器幅值更新
            end
            % aomig
            ao.para(2,2)=ao.para(2,2)+ao.para_dot(2,2)*dt;          %2阶振荡器频率更新
            if(ao.para(2,2)<0)
                ao.para(2,2)=0;             %频率非负限制
            end
            w_e=footp.freg(i)-ao.para(2,2);             %频率误差
            %w_e=0;
            ao.para(2,2)=ao.para(2,2)+5*w_e*dt;             %频率快速跟踪
            for j=3:order
            ao.para(j,2)=ao.para(2,2)*(j-1);                %高阶振荡器频率为基频整数倍
            end
            %fai
            for j=2:order
            ao.para(j,3)=ao.para(j,3)+ao.para_dot(j,3)*dt;          %j阶振荡器相位更新
            end    
            % obtain ao phase but starting by foot pressure
            
    else
            ao.y=0;
            ao.para_dot=zeros(order,3);         %振荡器参数变化率置零
    end
            
    % 
    gait.p_pure(i)=ao.para(2,3);                %步态相位纯值
    gait.p(i)=mod(ao.para(2,3),2*pi);           %步态相位mod 2pi
    gait.w(i)=ao.para(2,2);                 %步态频率
    gait.q(i)=ao.y;                         %估计关节角度
    gait.a0(i)=ao.para(1,1);                %振荡器0阶幅值
end


%% Predictive Locomotion Mode Recognition and Accurate Gait Phase Estimation for Hip Exoskeleton on Various Terrains
ep=0;               %相位误差初始化
phi=0;              %估计相位初始化
dphi=0;                 %相位变化率初始化
errorp=zeros(len,1);        %相位误差记录
delphi=zeros(len,1);            %估计相位记录
ddelphi=zeros(len,1);          %相位变化率记录
pcrt=zeros(len,1);             %相位校正值记录
kphi=3;                         %相位校正增益
count=0;                        %无步态触发计数器
h=0.3;                          %衰减系数
for i=1:len
    if(i>1)
        if(gflag(i)==0)             %步态停止时
            ep=0;    
        else
            if(footp.pha(i-1)==0&&footp.pha(i)==1)  %   步态触发时
                count=0;                        %无步态触发计数器清零
                ep=mod(gait.p_pure(i),2*pi);        %获取当前纯步态相位作为相位误差
                if(ep>pi)                   %相位误差调整到[-pi,pi]
                    ep=ep-pi;               %调整到[-pi,pi]
                end
                dphi=kphi*(ep-phi)*exp(0);          %相位变化率更新
            else
                count=count+1;                      %无步态触发计数器加一
                dphi=kphi*(ep-phi)*exp(-h*(count)*dt);              %相位变化率更新，指数衰减
                
            end
            phi=phi+dphi*dt;                %相位更新
            
        end
    end
    errorp(i)=ep;               %相位误差记录
    ddelphi(i)=dphi;                %相位变化率记录
    delphi(i)=phi;              %估计相位记录
    pcrt(i)=gait.p_pure(i)-delphi(i);               %相位校正值记录
end

gait.pcrt=pcrt;

figure
x1=subplot(3,1,1);
hold on
plot(t,qf_nb)
plot(t,gait.q)
plot(t,gflag)
plot(t,gait.a0)
% plot(t,rhip,t,lhip)
legend('q','q esti','flag','a0')
x2=subplot(3,1,2);
hold on
plot(t,footp.pha*2*pi)
plot(t,footp.freg)
plot(t,gait.w)
plot(t,gait.p)
legend('foot pha','foot freg','w','p')
x3=subplot(3,1,3);
plot(t,Euler)
linkaxes([x1,x2,x3],'x')



figure
x1=subplot(2,1,1);
plot(t,errorp)
hold on
plot(t,delphi)
plot(t,ddelphi/10)
legend('ep','phi','dphi')
x2=subplot(2,1,2);
plot(t,mod(pcrt,2*pi))
hold on
plot(t,footp.p)
%plot(t,gait.p)
legend('crt p','foot p')
linkaxes([x1,x2],'x')



end