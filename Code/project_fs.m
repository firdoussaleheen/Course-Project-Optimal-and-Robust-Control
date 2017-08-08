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
leadcomp = tf([-4.13 -22.07],[1 0]);
sys3 = series(leadcomp,engine);
sys4 = parallel(sys3,wind);
olsyscomp = series(sys4,aircraft);
figure(3);clf;
nyquist(olsyscomp)
figure(4);clf;
margin(olsyscomp);

clsyscomp = feedback(olsyscomp,1);

ltiview(clsyscomp,'r')
%step(clsyscomp,'r',clsyscomp2,'b--')