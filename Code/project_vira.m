clear all;
clc;
%%OLTF (K=2)
gaincomp = 2;
engine = tf(100,[1 10]);
wind = tf(-1,[1 0]);
aircraft = tf(40, [1 20 0]);
 
sys1 = series(gaincomp,engine);
sys2 = parallel(sys1, wind);
olsys = series(sys2, aircraft);
 	
figure(1);clf;
nyquist(olsys)
xlim([-10 30]);
%ylim([-40 40]);
grid on
%%
%find margin
figure(2);clf;
margin(olsys)
 %%
leadcomp = tf([1 15.77],[1 22.88]);
sys3 = series(leadcomp,engine);
sys4 = parallel(sys3,wind);
olsyscomp = series(sys4,aircraft);
figure(3);clf;
nyquist(olsyscomp)
figure(4);clf;
margin(olsyscomp);
%% 
leadcomp2 = tf([1 15.77],[1 10.27]);
sys5 = series(leadcomp,leadcomp2)
sys6 = series(sys5,engine);
sys7 = parallel(sys6,wind);
olsyscomp2 = series(sys7,aircraft)
figure(5);clf;
nyquist(olsyscomp2)
figure(6);clf;
margin(olsyscomp2);
 
%%
figure (7);clf;
bode(olsyscomp,'r')
hold on;
bode(olsyscomp2,'b--')
%%
clsyscomp = feedback(olsyscomp,1);
clsyscomp2 = feedback(olsyscomp2,1);
figure (8);clf;
%ltiview(clsyscomp,'r',clsyscomp2,'b--')
step(clsyscomp,'r',clsyscomp2,'b--')
%%
fb = bandwidth(clsyscomp)
fb2 = bandwidth(clsyscomp2)
