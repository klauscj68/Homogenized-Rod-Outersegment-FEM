close all


load SPR_mouse_wt
drop_mouse_wt=drop;
time_mouse=time_downsample;

load SPR_mouse_Caclamp
drop_mouse_Caclamp=drop;



load SPR_salamander_wt
drop_salamander_wt=drop;
time_salamander=time_downsample;


load SPR_salamander_Caclamp
drop_salamander_Caclamp=drop;




figure(1)
hold on
a=plot(time_mouse,drop_mouse_wt,'-b');
b=plot(time_mouse,drop_mouse_Caclamp,'-r');
set([a b],'linewidth',2);
aa=xlabel('time [s]');
bb=ylabel('drop [%]');
set([aa bb],'fontsize',20);
legh=legend([a b],'wt','Ca clamp');
set(legh,'fontsize',16);
c=title('Mouse');
set(c,'fontsize',25);

print -depsc mouse_SPR




figure(2)
hold on
a=plot(time_salamander,drop_salamander_wt,'-b');
b=plot(time_salamander,drop_salamander_Caclamp,'-r');
set([a b],'linewidth',2);
aa=xlabel('time [s]');
bb=ylabel('drop [%]');
set([aa bb],'fontsize',20);
legh=legend([a b],'wt','Ca clamp');
set(legh,'fontsize',16);
c=title('Salamander');
set(c,'fontsize',25);

print -depsc salamander_SPR

