clear 

root= "/Users/horay/CodesArchiv/rns_v2/Cont/";
dir = strcat( root, "J1.8_seq_MPA1.dat" );

A = importdata(dir);

e = A(:,1); e_ = [1.969141270E+15 2.179812359E+15];
rp= A(:,2); rp_= [9.404762707E-01 9.417196558E-01];
m = A(:,4); m_ = [2.510220273E+00 2.510247776E+00];
mb= A(:,3); mb_= [3.050007310E+00 3.050029154E+00];
Eb= mb-m;

fig=figure(2); hold on; box on;
%plot(e,Eb/max(Eb),'-','LineWidth',2);
plot(e,mb ,'.');
plot(e_,mb_ ,'.','MarkerSize',21);
%plot(e(end-2),m(end-2),'.','MarkerSize',20);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;

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

