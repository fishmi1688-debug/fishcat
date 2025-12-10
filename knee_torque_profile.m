function Knee=knee_torque_profile(gflag,footp,par,footF,kneepos,lknee,knee_vel)

dt=0.01;
len=length(gflag);
kneeAss=zeros(len,1);
kneePd=zeros(len,1);
t=0:dt:dt*(len-1);
lenp=length(par.t);

ass.fle_max_amp_=par.fle_max_amp_(1,2);
ass.fle_duration_=par.fle_duration_(1,2);
ass.ext_max_amp_=par.ext_max_amp_(1,2);
ass.ext_duration_=par.duration_(1,2)-par.fle_duration_(1,2);

%% 

kneeflag1=zeros(len,1);
kneeflag2=zeros(len,1);
kneeflag=zeros(len,1);
kneeflagf=zeros(len,1);


kneeBias=kneepos(1)-lknee(1);
knee_d=par.knee_pos_(1)+kneeBias;


for i=2:len
    
    if(gflag(i)==1)
        if(footp.pha(i-1)==0&&footp.pha(i)==1)
            for j=2:lenp
                if(t(i)>=par.t(j-1)&&t(i)<par.t(j))
                    % hip amplitute
                    knee_d=par.knee_pos_(j)+kneeBias;
                    break;
                end
            end
        end
        % obtain the foot phase
        p=footp.p(i);
        if(footp.pha(i-1)==0&&footp.pha(i)==1)
            % obtain the assistance parameter
            for j=2:lenp
                if(t(i)>=par.t(j-1)&&t(i)<par.t(j))
                    % hip amplitute
                    ass.fle_max_amp_=par.fle_max_amp_(j,2);
                    ass.fle_duration_=par.fle_duration_(j,2);
                    ass.ext_max_amp_=par.ext_max_amp_(j,2);
                    ass.ext_duration_=par.duration_(j,2)-par.fle_duration_(j,2);
                    break;
                end
            end
        end
        if(p<=ass.fle_duration_)
        p_=p/ass.fle_duration_;
        ass.tor=ass.fle_max_amp_*sin(p_*pi);
        kneeflag1(i)=1;
        pd_oloop=0;
        elseif(p>ass.fle_duration_&&p<(ass.fle_duration_+ass.ext_duration_))
         p_=(p-ass.fle_duration_)/ass.ext_duration_;
         ass.tor=ass.ext_max_amp_*sin(p_*pi);
         kneeflag1(i)=1;   
         %% knee position constraint
         kp=30*p_^2;
         kd=1*p_;
         pd_oloop=kp*(knee_d-kneepos(i))+kd*(knee_d-knee_vel(i));
        else
         ass.tor=0;
         kneeflag1(i)=0;
         pd_oloop=0;
        end
        
        %%
        if(footp.pha(i)==1)
            kneeflag2(i)=1;
        else
            kneeflag2(i)=0;
        end
        
    else
        ass.tor=0;  
        pd_oloop=0;
        kneeflag1(i)=0;
        kneeflag2(i)=0;
    end
    
    %% 
    kneeAss(i)=ass.tor;
    kneePd(i)=pd_oloop;
    %% 
    if(kneeflag1(i)==0||kneeflag2(i)==0)
    kneeflag(i)=0;
    else
    kneeflag(i)=1;   
    end
    % 0-1 position mode to force mode 
    if(kneeflag(i)>kneeflagf(i))
        kneeflagf(i)=kneeflagf(i-1)+0.5;
    %1-0 force mode to position mode    
    elseif(kneeflag(i)<kneeflagf(i))
        kneeflagf(i)=kneeflagf(i-1)-0.05;
    end
    %
    if(kneeflagf(i)>=1)
    kneeflagf(i)=1;
    elseif(kneeflagf(i)<=0)
    kneeflagf(i)=0;
    end
end

Knee.Ass=kneeAss;
Knee.flag=kneeflagf;

figure
x1=subplot(3,1,1);
plot(t,footp.p,'black','LineWidth',1)
legend('p','interpreter','latex')
ylabel('phase', 'interpreter','latex')
set(gca,'FontSize',14)
x2=subplot(3,1,2);
plot(t,kneeAss,'black','LineWidth',1)
hold on
plot(t,kneePd,'blue','LineWidth',1)
plot(t,kneeflagf,'red','LineWidth',1)
legend('knee torque','knee open','knee flag','interpreter','latex')
set(gca,'FontSize',14)
ylabel('Nm', 'interpreter','latex')

x3=subplot(3,1,3);
plot(t,footp.pha,'black',t,footF/600,'red','LineWidth',1)
legend('p','foot flag','interpreter','latex')
set(gca,'FontSize',14)
xlabel('time (s)', 'interpreter','latex')

linkaxes([x1,x2,x3],'x')


figure
hold on
plot(t,kneeAss)
plot(t,kneeflagf)
plot(t,kneePd)
plot(t,kneepos)
plot(t,knee_d*ones(len,1))
legend('ass','flag','open pd','pos','The')

end
