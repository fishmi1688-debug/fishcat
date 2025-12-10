function gait_phase_generation()
% 函数说明：步态相位生成主函数
clear all
%
load('xxx-itf4.mat')
% 数据提取
t=itf4.data.t;   %时间
lhip=itf4.data.sang(:,1);          %左髋关节角度
lknee=itf4.data.sang(:,2);         %左膝关节角度
rhip=itf4.data.sang(:,3);          %右髋关节角度
footF=itf4.data.footRT;            %脚底压力
Euler=itf4.data.Euler;             %欧拉角
kneepos=itf4.data.mpos(:,2);       %膝关节位置
par=itf4.PAR;                      %参数

knee_cur=itf4.data.mtor(:,2);      %膝关节电流
knee_vel=gradient(itf4.data.mpos(:,2))*100;    %膝关节速度

% bias estimation
preproc=bias_estimation(rhip,lhip);         %预处理
% start stop calssification
gflag=gait_flag_classification(t,preproc,footF);    %步态开始停止分类
% phase by foot pressure
footp=GaitPExtraction(t,footF,gflag);            %脚底压力相位提取
% real phase identification
gaitp_true=truePhaseObtainment(t,footp,gflag,footF);     %真实步态相位识别
% gait phase estimation
gaitp=adaptive_oscillators(preproc,gflag,footp,Euler);   %步态相位估计
% hip torque generation
hipAss=hip_torque_profile(gflag,footp,par);      %髋关节助力生成
% knee torque generation
kneeA=knee_torque_profile(gflag,footp,par,footF,kneepos,lknee,knee_vel);        %膝关节助力生成
% knee disturbance observer
dist=disturbance_observer(knee_cur,knee_vel);       %膝关节扰动观测
% knee joint command
knee_cmd=knee_joint_current_cmd(kneeA,kneepos,dist,knee_vel,par,lknee,footp);         %膝关节电流指令生成
% 



end

