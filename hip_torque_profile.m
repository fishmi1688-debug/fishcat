function hipAss=hip_torque_profile(gflag,footp,par)             %髋关节助力生成

dt=0.01;                        %时间步长
len=length(gflag);              %数据长度
hipAss=zeros(len,1);            %髋关节助力初始化
t=0:dt:dt*(len-1);              %时间序列
lenp=length(par.t);             %参数长度

ass.fle_max_amp_=par.fle_max_amp_(1,1);             %髋关节屈曲最大幅值
ass.fle_duration_=par.fle_duration_(1,1);           %髋关节屈曲持续时间
ass.ext_max_amp_=par.ext_max_amp_(1,1);             %髋关节伸展最大幅值
ass.ext_max_amp_=-4;                            %髋关节伸展最大幅值
ass.ext_duration_=2*pi-0.2*pi-ass.fle_duration_;            %髋关节伸展持续时间
for i=2:len                                     
    
    if(gflag(i)==1)
        p=footp.p(i);               %当前足相位
        if(footp.pha(i-1)==0&&footp.pha(i)==1)         %    步态触发时
            for j=2:lenp
                if(t(i)>=par.t(j-1)&&t(i)<par.t(j))                 
                    % hip amplitute
                    ass.fle_max_amp_=par.fle_max_amp_(j,1);         %髋关节屈曲最大幅值
                    ass.fle_duration_=par.fle_duration_(j,1);           %髋关节屈曲持续时间
                    ass.ext_max_amp_=par.ext_max_amp_(j,1);             %髋关节伸展最大幅值
                    ass.ext_max_amp_=-4;                            %髋关节伸展最大幅值
                    ass.ext_duration_=2*pi-0.1*pi-ass.fle_duration_;            %髋关节伸展持续时间
                    break;
                end
            end
        end
        if(p<=ass.fle_duration_)
        p_=p/ass.fle_duration_;                                     %归一化屈曲相位
        ass.tor=ass.fle_max_amp_*sin(p_*pi);                        %屈曲助力
        elseif(p>ass.fle_duration_&&p<(2*pi-0.1*pi))                %伸展助力
         p_=(p-ass.fle_duration_)/ass.ext_duration_;                %归一化伸展相位
         ass.tor=ass.ext_max_amp_*sin(p_*pi);                        %伸展助力
        else
         ass.tor=0;                                                 %无助力
        end
        
    else
        ass.tor=0;                                                   %无助力
    end
    hipAss(i)=ass.tor;
    
end


figure
x1=subplot(2,1,1);
plot(t,footp.p_c,'black','LineWidth',1)
legend('p','interpreter','latex')
ylabel('phase', 'interpreter','latex')
set(gca,'FontSize',14)

x2=subplot(2,1,2);
plot(t,hipAss,'red','LineWidth',1)
legend('hip torque','interpreter','latex')
xlabel('time (s)', 'interpreter','latex')
ylabel('Nm', 'interpreter','latex')
set(gca,'FontSize',14)
linkaxes([x1,x2],'x')


end
