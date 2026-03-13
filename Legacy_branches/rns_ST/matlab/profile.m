clear

dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.49658_B1.60E+01_mphi5.00E-02_rhoc1.365E+15_sphic0.869.dat";
dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.50048_B2.20E+01_mphi1.00E-01_rhoc1.540E+15_sphic0.962.dat";
dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.51006_B4.00E+01_mphi2.00E-01_rhoc1.730E+15_sphic1.067.dat";
dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.52049_B6.80E+01_mphi3.00E-01_rhoc1.782E+15_sphic1.105.dat";
%dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.52306_B1.40E+02_mphi5.00E-01_rhoc1.730E+15_sphic1.060.dat";
%dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.54283_B5.20E+02_mphi1.00E+00_rhoc1.662E+15_sphic1.050.dat";
%dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.54048_B4.00E+03_mphi3.00E+00_rhoc1.609E+15_sphic0.992.dat";
%dir = "/Users/horay/CodesArchiv/rns_ST/Cont/PS_J_0.00_Mb1.56426_B4.80E+04_mphi1.00E+01_rhoc1.597E+15_sphic1.031.dat";
dir = "/Users/horay/CodesArchiv/rns_ST/Res/res.dat";

A = importdata(dir," ",1);

res = sqrt( length(A.data(:,1)) );

SDIV = res;
MDIV = res;
s = A.data(1:MDIV:end-MDIV+1,1);
m = A.data(1:MDIV,2);

r = s./(1-s);
x = m.*r';
y = sin(acos(m)).*r';

cnt = 1;
for ss=1:SDIV
    for mm=1:MDIV
        alp(ss,mm)=A.data(cnt,3);
        gam(ss,mm)=A.data(cnt,4);
        rho(ss,mm)=A.data(cnt,5);
        ww (ss,mm)=A.data(cnt,6);
        prs(ss,mm)=log10( A.data(cnt,7));
        hh (ss,mm)=exp( A.data(cnt, 9) );
        rss(ss,mm)=log10( A.data(cnt,10));
        v2 (ss,mm)=sqrt( A.data(cnt,11) );
        Om (ss,mm)=A.data(cnt,12);
        sphi(ss,mm)=( A.data(cnt,13) );
        cnt = cnt+1;
    end
end

log_a = (rho+gam)/2;
fig=figure(103); fig.Position=[200 200 750 610];
p1 = [0.12,0.11,0.78,0.8]; subplot('Position',p1); hold on; box on;
phi=exp(sphi.^2/2);

%h=pcolor(y,x,hh'.*exp(log_a'-sphi'.^2/4)); 
h=pcolor(y,x,sphi'); tle = "$\varphi$";
set(h, 'EdgeColor', 'none'); 
cmap = load('/Users/horay/Data4Projects/gist_ncar_256.txt'); colormap(cmap);
%cmap = load('/Users/horay/Data4Projects/inferno_256.txt'); colormap(cmap);
colorbar; hcb=colorbar; set(hcb,'TickLabelInterpreter','latex');

L = 2.5; xlim([0 L]); ylim([0 L]);

xlabel("$x$ ($r_e$)",'Interpreter','latex');
ylabel("$z$ ($r_e$)",'Interpreter','latex');

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;

title(tle,'Interpreter','latex','FontSize',34);
%
set(fig,'Units','Inches'); pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig, '/Users/horay/Data4Projects/crazy/plots/tmp03','-dpdf','-r0');
%}
%%
clear
s = linspace(0,1-1e-10,241);
mphi = 0.1;

semilogy(s,s./(1-s),'-o');
%ylim([1 mphi*10*5])
%%
clear
mphi = 1e-1/1.4769994423016508;

figure(1); hold on; box on;
dir = strcat("/Users/horay/CodesArchiv/rns_ST/Debug/1d.dat");
A = importdata(dir);
s = A(:,1);
r = s./(1-s);
sphi = A(:,12);

plot(r, sphi ,'-o','LineWidth',2); xlim([0 .12]);
%plot(r, sphi , '-o','LineWidth',2); set(gca,'yscale','log'); %set(gca,'xscale','log'); xticks(logspace(-3,3,7))

xlabel("$s$",'Interpreter','latex');
ylabel("$\varphi$",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',24); ax = gca; ax.LabelFontSizeMultiplier = 1.2;