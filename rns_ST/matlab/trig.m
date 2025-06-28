clear 
clr = {[0.00,0.45,0.74],[0.85,0.33,0.10],[0.93,0.69,0.13],[0.49,0.18,0.56], ...
    [0.30,0.75,0.93],[0.78,0.54,0.86],[0.47,0.67,0.19,0.7],...
    [0.00,0.45,0.74],[0.85,0.33,0.10],[0.93,0.69,0.13],[0.49,0.18,0.56]};

root= "/Users/horay/CodesArchiv/rns_ST/trig/";

figure(38); hold on; box on; grid on;
for J = 7: 10
    dir = strcat( root, "mod_bessel_",num2str(J*2),".dat" );
    
    A = importdata(dir);
    
    r = A(:,1);
    k = A(:,2);
    i = A(:,3);
    
    %plot(r,i.*k,'color',clr{J+1},'LineWidth',2);
end
%plot(r, ( (r.^2+3)-3*r.*cosh(r)./sinh(r) ) ./ r.^3 )


s = linspace(1e-5,1,1200);
z = s./(1-s);
for nu = 7: 10
    I(nu+1,:) = besseli(nu*2+.5,z)*sqrt(pi/2)./sqrt(z);
    K(nu+1,:) = besselk(nu*2+.5,z)*sqrt(2/pi)./sqrt(z);
    %plot(z,K(nu+1,:),'-.','LineWidth',2);
    plot(z,I(nu+1,:),'o','color',clr{nu+1});
    %plot(z,I(nu+1,:).*z./sinh(z).*(K(nu+1,:).*z.*exp(z)),'o','color',clr{nu+1});
end

set(gca,'yscale','log'); 
set(gca,'xscale','log');
xlim([1e-2 1e6]);
%ylim([1e-7 1]);
