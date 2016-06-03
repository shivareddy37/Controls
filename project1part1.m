A= [-1 250 0 0 0 0 0;...
    -250 -1 0 0 0 0 0;...
    0 0 -1 190 0 0 0;...
    0 0 -190 -1 0 0 0;...
    0 0 0 0 -1 150 0;...
    0 0 0 0 -150 -1 0;...
    0 0 0 0 0 0 -300];
B= [21.1;...
    7.63;...
    7.65;...
    2.17;...
    6.45;...
    1.52;...
    25]; 
C= [7.63 21.1 2.17 7.65 1.52 6.45 25];
D= 0;
format short
n=size(A); 
olpoles=eig(A);
clpoles=[-50-250i; -50+250i; -38-190i; -38+190i; -30-150i; -30+150i; -300];
disp('the open loop pole loaction are:');
disp(olpoles);
disp('the desired closed loop pole loaction are:');
disp(clpoles);
olcharcoff=poly(olpoles);
clcharcoff=poly(clpoles);
disp('the coffecients of open loop characterstic equation are :');
disp(olcharcoff);
disp('the coffecients of desired closed loop characterstic equation are :');
disp(clcharcoff);

% the vale of kbar
kbar = clcharcoff(2:end) - olcharcoff(2:end);
disp('the kbar matrix is : ');
disp(kbar);

%evaluating the transformation matrix

% controllability matrix

Qc=zeros(n);
Qc(:,1)=B;
for k=2:n
    Qc(:,k)=A*Qc(:,k-1);
end
disp(' the controlibility matrix is :');
disp(Qc);
% the toeplitz matrix

col=zeros(n,1);
col(1) = 1;
Vinv = toeplitz( col, olcharcoff(1:end-1) );
disp('the teoplitz matrix is :');
disp(Vinv);

%the transformation matrix 

T=inv(Vinv)/Qc;
disp('the transformation matrix is :');
disp(T);
disp('THE REQUIRED GAIN MATRIX TO PUT THE OPEN LOOP POLES AT THE DESIRED CLOSED LOOP LOCATION IS :');
K=kbar*T;
disp(K);
norm_K=norm(K);
disp('the norm of the gain matrix using ackman"s formula is');
disp(norm_K);
%verification of eigon values of closed loop

Cleig=eig(A-B*K);
disp('the closed loop eigon values are ');
disp(Cleig);

%stumulation for closed loop homegeneous response


t = 0:0.001:5;
u = zeros(size(t));
x0=[-1.94;-1.65;-.78;.54;-.72;.31;-.87];

%Simulation
%Closed loop
figure('Name','Closed Loop homogeneous response ')
cl_system = ss(A-B*K,B,C,0);
lsim(cl_system,u,t,x0);
xlabel('Time (sec)')
ylabel('Amplitude')

%Open loop
figure('Name','Open Loop homogeneous response')
ol_system = ss(A,B,C,0);
lsim(ol_system,u,t,x0);
xlabel('Time (sec)')
ylabel('Amplitude')