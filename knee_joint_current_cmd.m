function knee_cmd=knee_joint_current_cmd(kneeA,kneepos,dist,knee_vel,par,lknee,footp)

lenp=length(par.t);

len=length(kneepos);
knee_loop=zeros(len,1);
kp_loop=zeros(len,1);
kd_loop=zeros(len,1);
pd_loop=zeros(len,1);
kp=30;
kd=2;
dt=0.01;
t=0:dt:dt*(len-1);
kneeBias=kneepos(1)-lknee(1);
knee_d=par.knee_pos_(1)+kneeBias;
for i=2:len
        if(footp.pha(i-1)==0&&footp.pha(i)==1)

            for j=2:lenp
                if(t(i)>=par.t(j-1)&&t(i)<par.t(j))
                    % hip amplitute
                    knee_d=par.knee_pos_(j)+kneeBias;
                    break;
                end
            end
        end
    kp_loop(i)=kp*(knee_d-kneepos(i));
    kd_loop(i)=kd*(knee_d-knee_vel(i));
    pd_loop(i)=kp_loop(i)+kd_loop(i);
    pd_loop(i)=(1-kneeA.flag(i))*pd_loop(i);
    knee_loop(i)=kneeA.flag(i)*kneeA.Ass(i)+(1-kneeA.flag(i))*(pd_loop(i)-dist(i));

end

knee_cmd=knee_loop;


figure
x1=subplot(3,1,1);
hold on
plot(t,kneeA.flag*3)
plot(t,kneeA.Ass)
plot(t,footp.p)
legend('kneeflag','ass','p')
x2=subplot(3,1,2);
hold 
plot(t,kp_loop)
plot(t,kd_loop)
plot(t,-dist)
plot(t,knee_loop)
legend('kp','kd','dist','loop')
x3=subplot(3,1,3);
hold on
plot(t,kneeA.flag)
plot(t,kneepos,t,knee_d*ones(len,1))
legend('p','pos a','pos d')
linkaxes([x1,x2,x3],'x')

end
