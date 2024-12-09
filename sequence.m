clear 

root= "/Users/horay/CodesArchiv/RNS_target/";
dir = strcat( root, "jbd1.dat" );

A = importdata(dir);

e = A(:,1);
rp= A(:,3);
ch= A(:,5);
m = A(:,6);
mb= A(:,7);

fig=figure(2); hold on; box on;
plot(e,m,'LineWidth',2);
%plot(e(end-2),m(end-2),'.','MarkerSize',20);


%%
clear

dir = "/Users/horay/CodesArchiv/RNS_target/check2D.dat";
A = importdata(dir);
SDIV = 1201;
MDIV = 901;


s = A(1:MDIV:end-MDIV,1);
m = A(1:MDIV,2);

r = s./(1-s);
x = m.*r';
y = sin(acos(m)).*r';

cnt = 1;
for ss=1:SDIV-1
    for mm=1:MDIV
        rho(ss,mm)=A(cnt,3);
        gam(ss,mm)=A(cnt,4);
        pre(ss,mm)=A(cnt,5);
        eng(ss,mm)=A(cnt,6);
        hh (ss,mm)=A(cnt,7);
        v2 (ss,mm)=A(cnt,8);
        ww (ss,mm)=A(cnt,10);
        cnt = cnt+1;
    end
end


fig=figure(102); fig.Position=[200 200 880 700]; box on;

h=pcolor(y,x,ww'); 
set(h, 'EdgeColor', 'none'); 
colormap hot; %clim([1 1.05]);
colorbar;
xlim([0 1.2]);
ylim([0 1.2]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;

