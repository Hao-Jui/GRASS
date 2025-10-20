clear

dir = "/Users/horay/CodesArchiv/rns_v2/eos/DD2_PT.dat";
A = importdata(dir," ",1);

e = A.data(:,1);
p = A.data(:,2);
h = log( A.data(:,3) );
n = A.data(:,4);

plot(p,e);

set(gca,'xscale','log')
set(gca,'yscale','log')