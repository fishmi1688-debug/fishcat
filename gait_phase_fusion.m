function gait_phase_fusion(preproc,gaitp,footp,gflag)

%%  Predictive Locomotion Mode Recognition and Accurate Gait Phase Estimation for Hip Exoskeleton on Various Terrains
len=length(gaitp.p);
errorp=zeros(len,1);
dt=0.01;
t=0:dt:dt*(len-1);
for i=1:len
    if(gflag==0)
        
    else
        % the fist pi phase
        if(footp.p(i)<pi)
        % error
        errorp(i)=footp.p(i)+pi-gaitp.p(i);
        errorp(i)=mod(errorp(i),2*pi);
        coefe=sin(footp.p(i));
        else
        errorp(i)=0;   
        coefe=0;  
        end
    end
end


figure
x1=subplot(3,1,1);
hold on
plot(t,preproc.qf_nb)
plot(t,gaitp.q)
legend('q','q esti')
x2=subplot(3,1,2);
hold on
plot(t,footp.p)
plot(t,footp.freg)
plot(t,gaitp.w)
plot(t,gaitp.p)
legend('foot p','foot freg','w','p')
x3=subplot(3,1,3);
plot(t,errorp)
linkaxes([x1,x2,x3],'x')


end
