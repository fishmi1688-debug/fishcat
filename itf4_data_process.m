
function itf4=itf4_data_process(filename)

%
clear all
filename='28';
fileID = fopen([filename,'.dat']);
Data = fread(fileID);
fclose(fileID);

gap=44;

head=Data(1:gap:end)*256^3+Data(2:gap:end)*256^2+Data(3:gap:end)*256+Data(4:gap:end);

t=Data(5:gap:end)+Data(6:gap:end)*256+Data(7:gap:end)*256^2+Data(8:gap:end)*256^3;

index_data=find(head==hex2dec('FFFDFFFD'));
index_para=find(head==hex2dec('FFFEFFFE'));

% note that the time may be delied
t_data=(t(index_data)-t(1))*0.0001;
t_para=(t(index_para)-t(1))*0.0001;

dt=0.01;
t_data=0:dt:(length(t_data)-1)*dt;
t_data=t_data';

para=zeros(length(t_para)*gap,1);
for i=1:length(index_para)
    para((i-1)*gap+1:i*gap)=Data((index_para(i)-1)*gap+1:index_para(i)*gap);
end

data=zeros(length(t_data)*gap,1);
for i=1:length(index_data)
    data((i-1)*gap+1:i*gap)=Data((index_data(i)-1)*gap+1:index_data(i)*gap);
end

%%  data  
sang=[];
for i=1:3
    tmp=data(9+2*(i-1):gap:end)+data(10+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/1000;
    sang=[sang; tmp(1:2:end)'];
end

mpos=[];
for i=1:2
    tmp=data(15+2*(i-1):gap:end)+data(16+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/1000;
    mpos=[mpos; tmp(1:2:end)'];
end

stor=[];
for i=1:2
    tmp=data(19+2*(i-1):gap:end)+data(20+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    stor=[stor; tmp(1:2:end)'];
end

mtor=[];
for i=1:2
    tmp=data(23+2*(i-1):gap:end)+data(24+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    mtor=[mtor; tmp(1:2:end)'];
end

mtord=[];
for i=1:2
    tmp=data(27+2*(i-1):gap:end)+data(28+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    mtord=[mtord; tmp(1:2:end)'];
end

knee_flag=data(31:gap:end)/100;

p=data(32:gap:end)/40;

% SYNC
sync=[];
tmp=data(33:gap:end)+data(34:gap:end)*256;
tmp=single(typecast(int32(tmp),'int16'));
sync=[sync; tmp(1:2:end)'];

% F total
footRT=[];
tmp=data(35:gap:end)+data(36:gap:end)*256;
tmp=single(typecast(int32(tmp),'int16'));
footRT=[footRT; tmp(1:2:end)'];

% q0
q0=[];
tmp=data(37:gap:end)+data(38:gap:end)*256;
tmp=single(typecast(int32(tmp),'int16'))/10000;
q0=[q0; tmp(1:2:end)'];

% q1
q1=[];
tmp=data(39:gap:end)+data(40:gap:end)*256;
tmp=single(typecast(int32(tmp),'int16'))/10000;
q1=[q1; tmp(1:2:end)'];

% q2
q2=[];
tmp=data(41:gap:end)+data(42:gap:end)*256;
tmp=single(typecast(int32(tmp),'int16'))/10000;
q2=[q2; tmp(1:2:end)'];

% q3
q3=[];
tmp=data(43:gap:end)+data(44:gap:end)*256;
tmp=single(typecast(int32(tmp),'int16'))/10000;
q3=[q3; tmp(1:2:end)'];

for i=1:length(q0)
    Quat(i)=quaternion([q0(i),q1(i),q2(i),q3(i)]);
end

Eluer=eulerd(Quat,'ZXZ','frame');

DATA.t=t_data;
DATA.sang=sang';
DATA.mpos=mpos';
DATA.stor=stor';
DATA.mtor=mtor';
DATA.mtord=mtord';
DATA.knee_flag=knee_flag;
DATA.p=p;
DATA.sync=sync';
DATA.footRT=footRT';
DATA.Euler=Eluer;
DATA.q=[q0',q1',q2',q3'];



%% parameter
ext_max_amp_=[];
for i=1:2
    tmp=para(9+2*(i-1):gap:end)+para(10+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    ext_max_amp_=[ext_max_amp_; tmp(1:2:end)'];
end

duration_=[];
for i=1:2
    tmp=para(13+2*(i-1):gap:end)+para(14+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    duration_=[duration_; tmp(1:2:end)'];
end

fle_max_amp_=[];
for i=1:2
    tmp=para(17+2*(i-1):gap:end)+para(18+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    fle_max_amp_=[fle_max_amp_; tmp(1:2:end)'];
end

fle_duration_=[];
for i=1:2
    tmp=para(21+2*(i-1):gap:end)+para(22+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    fle_duration_=[fle_duration_; tmp(1:2:end)'];
end

bias_=[];
for i=1:2
    tmp=para(25+2*(i-1):gap:end)+para(26+2*(i-1):gap:end)*256;
    tmp=single(typecast(int32(tmp),'int16'))/100;
    bias_=[bias_; tmp(1:2:end)'];
end

tmp=(para(29:gap:end)+para(30:gap:end)*256);
tmp=single(typecast(int32(tmp),'int16'))/100;
knee_pos_=tmp(1:2:end)';

tmp=(para(31:gap:end)+para(32:gap:end)*256);
tmp=single(typecast(int32(tmp),'int16'))/100;
ao_theh=tmp(1:2:end)';

name=[char(para(33:gap:end)),char(para(34:gap:end)),char(para(35:gap:end))];
preference=para(34:gap:end);

PAR.t=t_para;
PAR.fle_max_amp_=fle_max_amp_';
PAR.ext_max_amp_=ext_max_amp_';
PAR.fle_duration_=fle_duration_';
PAR.bias_=bias_';
PAR.knee_pos_=knee_pos_';
PAR.ao_theh=ao_theh';
PAR.name=name';
PAR.preference=preference';

itf4.data=DATA;
itf4.PAR=PAR;





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


save lsl_itf4.mat itf4

end
