sys_model=importdata('cbeam_3x2_sysmodel_Ts20kHz_194pole__struc.mat');
red_model=importdata('cbeam_3x2_reduced_model_40pole_Ts20kHz__struc.mat');
sys_FRF_data=importdata('cbeam_3x2_frfdata_Ts20kHz_hlfadv_20kPts__struc.mat');
% importing data from 194 pole system model

A_sys=sys_model.a;
B_sys=sys_model.b;
C_sys=sys_model.c;
D_sys=sys_model.d;
Ts_sys=sys_model.ts;
sys_ss=ss(A_sys,B_sys,C_sys,D_sys,Ts_sys);

% importing data from reduced 40 pole model of the system

A=red_model.a;
B=red_model.b;
C=red_model.c;
D=red_model.d;
Ts=red_model.ts;

controller_ss=ss(A,B,C,D,Ts);

% design of R and Q matrix for design of State feedback gain matrix K and
% error feedbak matrix L

W=zeros(40,1);


%modal_form=canon(red_ss,'modal');
for k=1:2:33
    p=eig(A(k:k+1,k:k+1));
    W(k)=abs((1/Ts)*(log(p(1))));
    W(k+1)=abs((1/Ts)*(log(p(2))));
end

for k=35:40
    p=eig(A(k,k));
W(k)=abs((1/Ts)*(log(p)));
end
H=abs((300*2*pi)./((W*1i)+((300*2*pi))));
Q=diag(H);

%loop to evaluate the correct value of multiplicity factor for R
for k=20:0.5:100
    R_K=k*eye(2);
    R_L=k*eye(3);
K=dlqr(A,B,Q,R_K);
m=dlqr(A.',C.',Q,R_L);
L=m.';
s=eig(A-B*K-L*C+L*D*K);
if(s<(1-1e-7))
    break 
end
end
disp('the optimal value of multiplicity factor of R which makes the system agressive yet assure the stand alone contyroller is stable is ')
disp (k)


% stability of standalone controller

A_oc=(A-(B*K)-(L*C)+(L*D*K));
A_oc_eig=eig(A_oc);

disp('the maximum absolute value of poles  of stand alone controller using computed K and L matrix ')
format long
max(abs(A_oc_eig))
format short





% crating and importing FRF data 
freq=sys_FRF_data.Frequency;
sys_frf=sys_FRF_data.ResponseData;
H_2=ss(A_oc,L ,K,0,Ts);
obser_based_cont_frf=freqresp(H_2,freq);

comp_ol_FRF=zeros(size(obser_based_cont_frf,1),size(sys_frf,2),numel(freq));
for k=1:numel(freq)
    comp_ol_FRF(:,:,k)=obser_based_cont_frf(:,:,k)*sys_frf(:,:,k);
end

% plot of char loci of compensated open loop system
figure(1)
char_loci=zeros(size (obser_based_cont_frf,1),numel(freq));
for k=1:numel(freq)
    char_loci(:,k)=eig(comp_ol_FRF(:,:,k));
end
char_loci_mag=20*log10(abs(char_loci));
char_loci_phase = radtodeg(angle(char_loci));
plot(char_loci_phase(1,:),char_loci_mag(1,:),'green--*');
hold on
plot(char_loci_phase(2,:),char_loci_mag(2,:),'blue--.');
hold on
ngrid
grid on
title('Nichols Chart')
xlabel('Phase(in degrees)')
ylabel('Magnitude(in DB)')
hold off

% compensated open loop model

comp_ol_A=[A-(B*K)-(L*C)+(L*D*K) L*C_sys;zeros(194,40) A_sys];
comp_ol_B=[L*D_sys;B_sys];
comp_ol_C=[K zeros(2,194)];
comp_ol_D=[0 0;0 0];
comp_ol_sys=ss(comp_ol_A,comp_ol_B,comp_ol_C,comp_ol_D,Ts);

% design of entire closed loop system

A_cl=[A_sys -B_sys*K;L*C_sys A_oc-L*D_sys*K];
B_cl=[B_sys;B-L*D-L*D_sys];
C_cl=[C_sys -D_sys*K];
D_cl=D_sys;

cl_sys=ss(A_cl,B_cl,C_cl,D_cl,Ts);
A_cl_eig=eig(A_cl);
disp('the maximum absolute value of the poles of this 234 closed loop model is');
format long
max(abs(A_cl_eig))
format short
cl_sys_frf=freqresp(cl_sys,freq);

% plot of MIMO root loci , compensated open loop and closed loop poles ,
% intended poles location
n = 20000;
ro=logspace(-3,6,n);
closed_loop_poles=zeros(234,n);

tic 
for i=1:n
    closed_loop_poles(:,i)=tzero(comp_ol_A,comp_ol_B,comp_ol_C,((1/ro(:,i))*eye(2,2))+(comp_ol_D));

disp(i)
end
toc

transmision_zeros=tzero(comp_ol_sys);
comp_ol_poles=eig(comp_ol_A);

figure(2)

plot(comp_ol_poles(:),'x r','markersize',6,'linewidth',3)
hold on
plot(transmision_zeros(:),'o g','markersize',10,'linewidth',2)
hold on
plot( closed_loop_poles(:),'^ c','markersize',6,'linewidth',3 )
hold on
plot(eig((A-B*K)),'s m','markersize',6,'linewidth',3 )
hold on
plot(eig((A-L*C)),'s m','markersize',8,'linewidth',2 )
zgrid_hires(2)
hold on
k=legend('open loop poles','transmission zeros','closed loop poles',' intended  cl poles by K ','intended  cl poles by L');
title('rootloci')
xlabel('real part')
ylabel('imaginery part')
hold off





% using matlab place function to genrate K and L matrix
Km=place(A,B,eig(A-B*K));
p_mat=place(A',C',eig(A-L*C));
Lm=p_mat.';
norm_Km=norm(Km);
norm_Lm=norm(Lm);
norm_K=norm(K);
norm_L=norm(L);
disp('the norm of Km and K matrix are as follows')
disp(norm_Km)
disp(norm_K)
disp('the norm of Lm and L matrix are as follows')
disp(norm_Lm)
disp(norm_L)

o=eig(A-B*Km-Lm*C+Lm*D*Km);

disp('the maximum absolute value of poles of standalone controller using matrix genrated by  matlab place function ')
disp(max(abs(o)))


% frf magnitude computation

uncomp_sys_frf=freqresp(controller_ss,freq);
mag_uncomp_sys=20*log10(abs(uncomp_sys_frf));
mag_cl_sys=20*log10(abs(cl_sys_frf));

% bode FRF magnitube  plots
figure(3)
% plot of uncompensated sys bode FRF magnitude
semilogx(freq,squeeze(mag_uncomp_sys(1,1,:)),' r')
hold on
semilogx(freq,squeeze(mag_uncomp_sys(2,1,:)),'g')
hold on
semilogx(freq,squeeze(mag_uncomp_sys(3,1,:)),'b')
hold on
semilogx(freq,squeeze(mag_uncomp_sys(1,2,:)),'c')
hold on
semilogx(freq,squeeze(mag_uncomp_sys(2,2,:)),'m')
hold on
semilogx(freq,squeeze(mag_uncomp_sys(3,2,:)),'y')
hold on

%plot of closed loop sys bode FRF magnitude
semilogx(freq,squeeze((mag_cl_sys(1,1,:))),'-*r')
hold on
semilogx(freq,squeeze(mag_cl_sys(2,1,:)),'-*g')
hold on
semilogx(freq,squeeze(mag_cl_sys(3,1,:)),'-*b')
hold on
semilogx(freq,squeeze(mag_cl_sys(1,2,:)),'-*c')
hold on
semilogx(freq,squeeze(mag_cl_sys(2,2,:)),'-*m')
hold on
semilogx(freq,squeeze(mag_cl_sys(3,2,:)),'-*y')
hold on

grid on
p=legend('uncomp(in1,out1)','uncomp(in1,out2)','uncomp(in1,out3)','uncomp(in2,out1)','uncomp(in2,out2)','uncomp(in2,out3)','cl(in1,out1)','cl(in1,out2)','cl(in1,out3)','cl(in2,out1)','cl(in2,out2)','cl(in2,out3)');
set(p,'location','northwest')
title('the bode frf magnitude odf closed loop system and uncompensated open loop system')

hold off





