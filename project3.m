%importing the data for the flexible system

sys=importdata('ee661_proj3_ss_model.mat');
A=sys.a;
B=sys.b;
C=sys.c;
D=sys.d;
n=size(A);

% analysing the open loop poles of the system
disp('the open loop poles of the system are at')
ol_poles=eig(A);
disp(ol_poles)

% the desired closed loop pole location

des_cl_poles=1.005.*[-65+5900i...
                 -65-5900i...
                 -3+180i...
                 -3-180i...
                 -4+130i...
                 -4-130i...
                 -5+50i...
                 -5-50i];
     
% design of state-feedbak matrix K
psi=zeros(11,8);
 for i=1:n
    p=( cat(2,( des_cl_poles(i)*eye(n)-A ),B));
    [U,S,V]=svd(p);
    psi(:,i)=V(:,9)+V(:,10)+V(:,11);

 end


K=psi(9:end,:)/( psi(1:8,:) );


% tests to check the correctness of gain matrix
norm_K=norm(K);
actual_cl_poles=eig(A-B*K);
disp('the closed loop poles as computed using the designed gain matrix')
disp(actual_cl_poles);

% desired poles for the observer-error feedback matrix L

design_poles_L=0.995.*[-65+5900i...
                       -65-5900i...
                       -3+180i...
                       -3-180i...
                       -4+130i...
                       -4-130i...
                       -5+50i...
                       -5-50i...
                       ];
                   
% design of observer error matrix L

psi1=zeros(11,8);
 for i=1:n
    p=( cat(2,( design_poles_L(i)*eye(n)-A' ),C'));
    [U,S,V]=svd(p);
    psi1(:,i)=V(:,9)+V(:,10)+V(:,11);

 end

p=psi1(9:end,:)/( psi1(1:8,:) );
L=p.';


% tests to check correctness observer error matrix L
norm_L=norm(L);
actual_poles_L=eig(A-L*C);
disp('the poles which we get using the computed observer error matrix')
disp(actual_poles_L)

%the stand alone contoller system
A_oc=A-B*K-L*C+L*D*K;
cont_ol_poles=eig(A_oc);
disp('open loop poles of the stand alone system')
disp(cont_ol_poles)

%  FRF of compensated open loop system

w=(1:100000);
G=ss(A,B,C,D);
H=ss(A_oc,L,K,0);
sys_FRF=freqresp(G,w);
cont_FRF=freqresp(H,w);

comp_ol_FRF=zeros(size(sys_FRF,1),size(cont_FRF,2),numel(w));
for k=1:numel(w)
    comp_ol_FRF(:,:,k)=cont_FRF(:,:,k)*sys_FRF(:,:,k);
end

% plot of char loci of compensated open loop system
figure(1)
char_loci=zeros(size (cont_FRF,1),numel(w));
for k=1:numel(w)
    char_loci(:,k)=eig(comp_ol_FRF(:,:,k));
end
char_loci_mag=20*log10(abs(char_loci));
char_loci_phase = radtodeg(angle(char_loci));
plot(char_loci_phase(1,:),char_loci_mag(1,:),'green--*');
hold on
plot(char_loci_phase(2,:),char_loci_mag(2,:),'blue--.');
hold on
plot(char_loci_phase(3,:),char_loci_mag(3,:),'red--d');
hold on

ngrid
grid on
title('Nichols Chart')
xlabel('Phase(in degrees)')
ylabel('Magnitude(in DB)')
hold off


% the closed loop system

A_cl=[A -B*K;L*C A-B*K-L*C];
B_cl=[B;B];
C_cl=[C -D*K];
D_cl=0;
cl_sys=ss(A_cl,B_cl,C_cl,D_cl);
cl_poles=eig(A_cl);
disp('the closed loop poles of 16 pole model')
disp(cl_poles);


% plot of both compensated open loop and closed loop system 

%plotting FRF magnitude of uncomepnsated open loop system
figure(2)

uncomp_ol_sys=zeros(size (cont_FRF,1),numel(w));
for k=1:numel(w)
    uncomp_ol_sys(:,k)=eig(sys_FRF(:,:,k));
end
uncomp_sys_FRF_mag=20*log10(abs(uncomp_ol_sys));


plot(w,uncomp_sys_FRF_mag(1,:),'green--*');
hold on
plot(w,uncomp_sys_FRF_mag(2,:),'blue--.');
hold on
plot(w,uncomp_sys_FRF_mag(3,:),'red--d');
hold on

%plotting entire closed loop system

G_cl_sys=cl_sys;
comp_cl_sys_FRF=freqresp(G_cl_sys,w);
cont_FRF=freqresp(H,w);



complete_cl_sys=zeros(size (comp_cl_FRF,1),numel(w));
for k=1:numel(w)
    complete_cl_sys(:,k)=eig(comp_cl_sys_FRF(:,:,k));
end
cl_sys_FRF_mag=20*log10(abs(complete_cl_sys));


plot(w,cl_sys_FRF_mag(1,:),'c--*');
hold on
plot(w,cl_sys_FRF_mag(2,:),'m--.');
hold on
plot(w,cl_sys_FRF_mag(3,:),'y--d');
hold on
set(gca,'xscale','log')
grid on
title('bode plot of open loop compensated system and entire closed loop system')
xlabel('frequency (in rad)')
ylabel('Magnitude(in DB)')
s=legend('direct path 1 uncomp sys','direct path 2 uncomp sys','direct path 3 uncomp sys','direct path 1 cl sys','direct path 2 cl sys','direct path 3 cl sys');
hold off

% grammian to test contolibility and observerablity of closed loop system

Wc=gram(cl_sys,'c');
Wo=gram(cl_sys,'o');

eig_controllable=eig(Wc);
eig_observerable=eig(Wo);
disp('eigen values of grammian controllable matrix')
disp(eig_controllable)
disp('eigen values of grammian observerable matrix')
disp(eig_observerable)




% use of place function

Km=place(A,B,des_cl_poles);
p_mat=place(A',C',design_poles_L);
Lm=p_mat.';
norm_K_m=norm(Km);
norm_L_m=norm(Lm);

% comaparison of norms

disp('norm of computed K matrix')
disp(norm_K)
disp('norm of K matrix from place function')
disp(norm_K_m)
disp('norm of computed L matrix')
disp(norm_L)
disp('norm of L matrix from place fucntion')
disp(norm_L_m)


% checking the stability of standalone controller with matlab genrated gain matrix

A_standalone_cont=A-B*Km-Lm*C+Lm*D*Km;
disp( 'the poles of the standalone controller using matlab gentared gain matrix')
m=eig(A_standalone_cont);
disp(m)




             
                
                