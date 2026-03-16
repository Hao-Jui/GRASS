clear

dir = "/Users/horay/CodesArchiv/rns_universal/eos/muses.dat";
A = importdata(dir," ",1);

e = A.data(:,1);
p = A.data(:,2);
h = log( A.data(:,3) );
n = A.data(:,4);

plot(p,e);

set(gca,'xscale','log')
set(gca,'yscale','log')