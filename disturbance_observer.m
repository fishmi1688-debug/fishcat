function dist=disturbance_observer(knee_cur,knee_vel)


knee_B=0.4;
knee_J=0.01;
knee_l=30.0;
knee_c=knee_l*knee_J;

sys=tf(1,[1/60 1]);
coef=c2d(sys,0.01,'tutin');
qfil.num=coef.Numerator{1};
qfil.den=coef.Denominator{1};

len=length(knee_cur);
dist=zeros(len,1);
dt=0.01;
t=0:dt:dt*(len-1);

for i=1:len
    vel=knee_vel(i);
    cur=knee_cur(i);
	in=(knee_B-knee_c)*vel-cur;
	in_f=lowpass_filter_current_knee(in,qfil);
	out=in_f+knee_c*vel;
    
    dist(i)=out;
end


figure
plot(t,-dist)
hold on
plot(t,knee_cur)
legend('dist','cur')
end



function y=lowpass_filter_current_knee(u,qfil)
persistent u_last y_last
if isempty(u_last)
    u_last=0;
end
if isempty(y_last)
    y_last=0;
end

% y=0.9802*y_last+0.009901*u_last+0.009901*u;
y=-qfil.den(2)*y_last+qfil.num(2)*u_last+qfil.num(1)*u;
u_last=u;
y_last=y;
end
