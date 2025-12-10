function data_visualization()

clear all
load('lsl_itf4.mat')

figure
x1=subplot(2,1,1);
plot(itf4.data.t,itf4.data.Euler(:,2:3))
legend('plantar flexion','rotation','interpreter','latex')
x2=subplot(2,1,2);
plot(itf4.data.t,itf4.data.footRT,'red',itf4.data.t, itf4.data.knee_flag*1000,'m')
hold on
legend('FT','Knee flag','interpreter','latex')
linkaxes([x1,x2],'x')

figure
x1=subplot(3,1,1);
plot(itf4.data.t, itf4.data.sang(:,1),'blue',itf4.data.t, itf4.data.sang(:,3),'black',itf4.data.t, itf4.data.p/2/pi,'blue',itf4.data.t, itf4.data.knee_flag,'m')
legend('right hip','left hip', 'p','knee flag','interpreter','latex')
title('ITF4','interpreter','latex')
x2=subplot(3,1,2);
plot(itf4.data.t,itf4.data.sang(:,3)-itf4.data.sang(:,1),'blue',itf4.data.t,itf4.data.sang(:,2),'black',itf4.data.t,itf4.data.knee_flag,'red')
legend('q','kneepos','kneeflag')
x3=subplot(3,1,3);
plot(itf4.data.t,itf4.data.mtord(:,1),'red',itf4.data.t,itf4.data.mtord(:,2),'--red',itf4.data.t,itf4.data.mtor(:,1),'black',itf4.data.t,itf4.data.mtor(:,2),'--black')
legend('hip tord','knee tord','hip tor','knee tor','interpreter','latex')
linkaxes([x1,x2,x3],'x')



%% torque sensor
figure
x1=subplot(2,1,1);
plot(itf4.data.t,itf4.data.stor(:,1),'red',itf4.data.t,itf4.data.mtor(:,1),'blue')
legend('hip tors','hip tor','interpreter','latex')
set(gca,'xtick',[])
ylabel('Nm')
set(gca,'fontsize',20)
title('Torque Sensor and Motor Current','interpreter','latex')
x2=subplot(2,1,2);
plot(itf4.data.t,itf4.data.stor(:,2),'red',itf4.data.t,itf4.data.mtor(:,2),'blue')
legend('knee tors','knee tor','interpreter','latex')
xlabel('time/s')
ylabel('N')
set(gca,'fontsize',20)
linkaxes([x1,x2],'x')


%% Hip angle, foot force, interactive force relationship
figure
x1=subplot(3,1,1);
plot(itf4.data.t, itf4.data.sang(:,1),'blue',itf4.data.t, itf4.data.sang(:,3),'black',itf4.data.t, itf4.data.p/2/pi,'blue',itf4.data.t, itf4.data.knee_flag,'m')
legend('right hip','left hip', 'p','knee flag','interpreter','latex')
title('ITF4','interpreter','latex')
x2=subplot(3,1,2);
plot(itf4.data.t,itf4.data.footRT,'red',itf4.data.t, itf4.data.knee_flag*1000,'m')
hold on
legend('FT','Knee flag','interpreter','latex')
x3=subplot(3,1,3);
plot(itf4.data.t,itf4.data.mtord(:,1),'red',itf4.data.t,itf4.data.mtord(:,2),'black')
legend('hip tord','knee tord','interpreter','latex')
linkaxes([x1,x2,x3],'x')




end