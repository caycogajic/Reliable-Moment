% Toy model
% Plot Figure 1

clear all; clc

Z = @(a,b) 1+3*exp(-a)+3*exp(-2*a+b)+exp(-3*a+3*b);
P0 = @(a,b) 1/Z(a,b);
P1 = @(a,b) exp(-a)/Z(a,b);
P2 = @(a,b) exp(-2*a+b)/Z(a,b);
P3 = @(a,b) exp(-3*a+3*b)/Z(a,b);
P = @(a,b) [P0(a,b),P1(a,b),P1(a,b),P1(a,b),P2(a,b),P2(a,b),P2(a,b),P3(a,b)];

P0_ri = @(a,b) 1/Z(a,b);
P1_ri = @(a,b) 1/Z(a,b);
P2_ri = @(a,b) 1/Z(a,b);
P3_ri = @(a,b) exp(.6)/Z(a,b);
P_ri = @(a,b) [P0_ri(a,b),P1_ri(a,b),P1_ri(a,b),P1_ri(a,b),P2_ri(a,b),P2_ri(a,b),P2_ri(a,b),P3_ri(a,b)];

P_rm = @(a,b) [P0(a,b),P1(a,b),P1(a,b),P1(a,b),P2(a,b),P2(a,b),P2(a,b),P3(a,b)];

figure, bar([P(1,1.2);P_ri(1,1.2);P_rm(1,1.2)]')
set(gca,'XTicklabel',{'000','001','010','100','110','101','011','111'})
set(gca,'FontSize',15)
xlabel('Spiking pattern')
ylabel('Probability')