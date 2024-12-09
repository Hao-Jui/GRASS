clear

l_uni      = 1.4769994423016508;

dir = "/Users/horay/Data4Projects/RNS_target/data/checkTOV.dat";
A    = importdata(dir);
r_is = A(:,1); r_e = max(r_is);
r    = A(:,2);
m    = A(:,3);
e    = A(:,4);
p    = A(:,5);
rho  = A(:,6)*1.6749286e-24;
nu   = A(:,7);
lam  = A(:,8);

fig=figure(102); fig.Position=[200 200 1050 700]; hold on; box on; grid on;
plot(r,exp(nu*2),'LineWidth',2);

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;

%%
clear

dir = "/Users/horay/CodesArchiv/rns_v2/Cont/J1.80_Mb1.80.dat";
%dir = "/Users/horay/CodesArchiv/rns_v2/Cont/check2D.dat";
A = importdata(dir);


res = 200;


SDIV = res+1;
MDIV = res+1;
s = A(1:MDIV:end-MDIV+1,1);
m = A(1:MDIV,2);

r = s./(1-s);
x = m.*r';
y = sin(acos(m)).*r';

cnt = 1;
for ss=1:SDIV
    for mm=1:MDIV
        alp(ss,mm)=A(cnt,3);
        gam(ss,mm)=A(cnt,4);
        rho(ss,mm)=A(cnt,5);
        ww (ss,mm)=A(cnt,6);
        prs(ss,mm)=log10( A(cnt,7));
        hh (ss,mm)=A(cnt, 9);
        rss(ss,mm)=log10( A(cnt,10));
        v2 (ss,mm)=sqrt( A(cnt,11) );
        %Om (ss,mm)=A(cnt,12);
        cnt = cnt+1;
    end
end
log_a = (rho+gam)/2;
fig=figure(102); fig.Position=[200 200 880 700]; box on;

h=pcolor(y,x,hh'); 
set(h, 'EdgeColor', 'none'); colormap jet; %clim([0 4]);
colorbar; hcb=colorbar; set(hcb,'TickLabelInterpreter','latex');
L = 1.2; xlim([0 L]); ylim([0 L]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;
%%
clear
dir = "/Users/horay/CodesArchiv/rns_v2/Omega.dat";
A = importdata(dir);
r = A(:,1); Omg = A(:,2); ww = A(:,3); wO = A(:,4);

plot(r,Omg/Omg(1),'.','MarkerSize',8); xlim([0 r(end)/2]);
%set(gca,'xscale','log');set(gca,'yscale','log');
%%
clear

dir = "/Users/horay/CodesArchiv/rns_v2/chech_uryu.dat";
A = importdata(dir);

dir = "/Users/horay/CodesArchiv/rns_v2/check_uryu_1D.dat";
B = importdata(dir);
SDIV = 501;
MDIV = 301;


s = A(1:MDIV:end-MDIV+1,1);
m = A(1:MDIV,2);

r = s./(1-s);
x = m.*r';
y = sin(acos(m)).*r';

cnt = 1;
for ss=1:SDIV
    for mm=1:MDIV
        Fj(ss,mm)=A(cnt,3);
        hh (ss,mm)=A(cnt,4);
        omg(ss,mm)=A(cnt,5);
        v2 (ss,mm)=A(cnt,6);
        ww (ss,mm)=A(cnt,7);
        cnt = cnt+1;
    end
end
fig=figure(102); fig.Position=[200 200 880 700]; box on;

h=pcolor(y,x,hh'); 
set(h, 'EdgeColor', 'none'); colormap jet; %clim([0.25 .3]);
colorbar; hcb=colorbar; set(hcb,'TickLabelInterpreter','latex');
L = 1.2; xlim([0 L]); ylim([0 L]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;
fig=figure(103); plot(B(:,1)./(1-B(:,1)),B(:,2)); xlim([0 L]);

%%
clear
dir = "/Users/horay/CodesArchiv/rns_v2/Cont/";
cnt = 1;
fold(cnt)="J1.8"; lbd(cnt)="J=1.8"; cnt=cnt+1;
%fold(cnt)="J1.5"; lbd(cnt)="J=1.5"; cnt=cnt+1;
%fold(cnt)="J0.8"; lbd(cnt)="J=0.8"; cnt=cnt+1;


fig=figure(102); fig.Position=[200 200 680 400]; box on; hold on;
for i=1:cnt-1
    A = importdata(strcat(dir,fold(i),"_seq_MPA1.dat") );
    ec= A(:,1);
    rep=A(:,2);
    Mb= A(:,3);
    M = A(:,4);
    re= A(:,5);
    Omg_c=A(:,6);
    plot(ec,Mb,'.','DisplayName',lbd(i),'LineWidth',2);
    plot([7.02088248066766625E+14 2.455211871E+15],[1.800005437 3.039497279],'.','MarkerSize',30,'HandleVisibility','off');
end
ylim([1.75 3.2]);
xlabel("$e$ (cgs)",'Interpreter','latex'); ylabel("$M_0$ ($M_\odot$)",'Interpreter','latex');
legend('Interpreter','latex','Location','southeast','box','on');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',22); ax = gca; ax.LabelFontSizeMultiplier = 1.2;
