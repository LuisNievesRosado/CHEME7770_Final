
eta = linspace(0,1,20);

m = length(eta);
OUT = zeros(1,m);

for k = 1:m
    k=k
OUT05(k) = MC(0.5,eta(k));
OUT1(k) = MC(1,eta(k));
OUT5(k) = MC(5,eta(k));
end
t = 18;
plot(eta(1:t),OUT05(1:t),eta(1:t-1),OUT1(1:t-1),eta(1:t),OUT5(1:t),'LineWidth',2)
ylim([-200,10])
legend('\kappa_s = 0.5','\kappa_s = 1','\kappa_s = 5')
xlabel('\eta')
ylabel('\nu_{actin}')