function gaitp_true=truePhaseObtainment(t,footp,gflag,footF)


len=length(t);          %数据长度
p=zeros(len,1);         %步态相位初始化

index=[];           %存储每个步态周期起始点索引

for i=1:len
    if(gflag(i)==0)         %非步态阶段
        
    else
        if(footp.p(i)==0)       %步态阶段且步态相位为0，记录索引
           index=[index,i];             %记录索引
        end     
    end
end


for i=1:length(index)-1             %每个步态周期内步态相位赋值
    cnt=index(i+1)-index(i);            %当前步态周期长度
    p(index(i):index(i+1)-1)=0:2*pi/cnt:(2*pi-2*pi/cnt);            %步态相位赋值
    
end


for i=1:len
    if(gflag(i)==0)         %非步态阶段，步态相位置0
     p(i)=0;                %步态相位置0
    end
end

gaitp_true.p=p;             %输出真实步态相位

p_d=zeros(len,1);           %基于脚底压力的步态相位修正初始化

for i=1:len
    if(gflag(i)==0)         %非步态阶段
     p_d(i)=0;              %步态相位置0
    else
        if(footp.p(i)>0.8*2*pi&&footp.coef(i)>0.1)      %步态阶段且脚底压力相位大于0.8*2*pi且脚底压力系数大于0.1
            p_d(i)=(footp.coef(i))*2*pi;            %步态相位置修正为脚底压力系数对应的步态相位
        else
            p_d(i)=0;           %否则步态相位置0
        end
    end
end




figure
x1=subplot(2,1,1);
hold on
plot(gaitp_true.p)
plot(footp.p)
plot(footp.coef*2*pi)
plot(footF/300)
legend('p true','p','coef','footF')
x2=subplot(2,1,2);
hold on
plot(gaitp_true.p)
plot(p_d)
legend('p true','pd')
linkaxes([x1,x2],'x')

end
