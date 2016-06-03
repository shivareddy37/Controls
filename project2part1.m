cbeam_plant=importdata('cbeam_3x2_sysmodel_Ts20kHz_194pole__struc.mat');
cbeam_controller=importdata('cbeam_controller_2x3_Ts20kHz_40pole__struc.mat');

%extracting the controller mtarix's

A_bar=cbeam_controller.a;
B_bar=cbeam_controller.b;
C_bar=cbeam_controller.c;
D_bar=cbeam_controller.d;
ts_bar=cbeam_controller.ts;

% extracting the plant matrix's
A=cbeam_plant.a;
B=cbeam_plant.b;
C=cbeam_plant.c;
D=cbeam_plant.d;
ts=cbeam_plant.ts;

% compensated open loop model
comp_ol_A=[A_bar B_bar*C;zeros(194,40) A];
comp_ol_B=[B_bar*D;B];
comp_ol_C=[C_bar D_bar*C];
comp_ol_D= D_bar*D;

comp_ol=ss(comp_ol_A,comp_ol_B,comp_ol_C,comp_ol_D,ts);

A_delta=comp_ol_A;
B_delta=comp_ol_B;
C_delta=comp_ol_C;
D_delta=comp_ol_D;
%GEVP
G=[comp_ol_A,comp_ol_B;comp_ol_C,comp_ol_D];
M=[eye(234,234) zeros(234,2);zeros(2,234) zeros(2,2)];
[evecs,evals]=eig(G,M);

transmission_zeros=diag(evals);
matlab_trans_zeros=tzero(comp_ol);

plot(transmission_zeros(:),'s r','markersize',6,'linewidth',3)
hold on
plot(matlab_trans_zeros(:) ,'^ g','markersize',6,'linewidth',3)
hold on
title('the plot of evluated and matlab transmission zeros')
s=legend('computed transmission zeroes','matlab transmission zeroes');
set(s,'location','southwest')
hold off

%closed loop at diffrent gain factor ro
n = 20000;
ro=logspace(-3,6,n);
closed_loop_poles=zeros(234,n);

tic 
for i=1:n
    closed_loop_poles(:,i)=tzero(A_delta,B_delta,C_delta,((1/ro(:,i))*eye(2,2))+(D_delta));

disp(i)
end
toc
cl_poles_unitygain = tzero( A_delta,B_delta,C_delta,eye(2,2)+(D_delta));
open_loop_poles=eig(comp_ol_A);

%plot of the root loci

figure(2)
plot(open_loop_poles(:),'x g','markersize',6,'linewidth',3)
hold on
plot(transmission_zeros(:),'o r','markersize',6,'linewidth',3)
hold on
plot( closed_loop_poles(:),'.','markersize',6,'linewidth',3 )
hold on
plot(cl_poles_unitygain(:),'s c','markersize',8,'linewidth',2 )
zgrid_hires(2)
hold on
k=legend('open loop poles','transmission zeros','closed loop poles','closed loop poles with unity gain');
set(k,'location','southwest')
title('rootloci')
xlabel('real part')
ylabel('imaginery part')
hold off
