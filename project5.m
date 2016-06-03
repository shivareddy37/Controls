
% importing data

sys=importdata('proj5_abcd_ts.mat');
input=importdata('proj5_input_sequence.mat');
noise=importdata('proj5_noise_sequences.mat');

A=sys.a;
B=sys.b;
C=sys.c;
D=sys.d;
ts=sys.ts;
V=noise.V;
W=noise.W;
Xi=noise.Xi;

R=0.003*eye(3);
Q=0.0025*eye(3);
sigma=0.00001*eye(10);


% initsilisation of system

X=randn(10,1);
Y=zeros(3,20000);
Z=zeros(3,20000);
P=zeros(10,10,20000);
X_cap=zeros(10,20000);
Y_cap=zeros(3,20000);
K_bar_k=zeros(10,3,20000);
norm_P=zeros(20000,1);
norm_K=zeros(20000,1);

P_init=100*eye(10,10);

K_0=(A*P_init*C')*((C*P_init*C')+inv(R));
P(:,:,1)=(A-(K_0*C))*P_init*(A-(K_0*C))'+(K_0*R*K_0')+(B*Q*B')+sigma;
norm_P(1,1) = norm(P(:,:,1));

%stimulation of the system

for k=1:1:19999
    X(:,(k+1))=A*X(:,k)+B*(input(:,k+1)+W(:,k+1))+Xi(:,k+1);
    Y(:,k)=C*X(:,k);
    Z(:,k)=C*X(:,k)+V(:,k+1);
    K_bar_k(:,:,k)=A*P(:,:,k)*C'/(C*P(:,:,k)*C'+R);
    P(:,:,k+1)=(A-K_bar_k(:,:,k)*C)*P(:,:,k)*(A-K_bar_k(:,:,k)*C)'+(K_bar_k(:,:,k)*R*(K_bar_k(:,:,k))'+B*Q*B'+sigma);
    X_cap(:,k+1)=(A-K_bar_k(:,:,k)*C)*X_cap(:,k)+B*input(:,k+1)+K_bar_k(:,:,k)*Z(:,k);
    Y_cap(:,k)=C*X_cap(:,k);
    norm_P(k+1,1)= norm(P(:,:,k+1));
    norm_K(k,1) = norm(K_bar_k(:,:,k));
end 
    
    


% plot

% true input and estimated input plots

figure,
plot(X(1,:),'k');
hold on
plot(X_cap(1,:),'b');
hold on
plot(X(2,:),'g');
hold on
plot(X_cap(2,:),'r');
hold on
xlabel('Time');
ylabel('State Vector');
legend('true state 1','Estimated state 1','true state 2','Estimated state 2');
title('State #1 and State #2 Vs Time plot');

figure,
plot(X(3,:),'k');
hold on
plot(X_cap(3,:),'b');
hold on
plot(X(4,:),'g');
hold on
plot(X_cap(4,:),'r');
hold on
xlabel('Time');
ylabel('State Vector');
legend('true state 3','Estimated state 3','true state 4','Estimated state 4');
title('State #3 and State #4 Vs Time plot');

figure,
plot(X(5,:),'k');
hold on
plot(X_cap(5,:),'b');
hold on
plot(X(6,:),'g');
hold on
plot(X_cap(6,:),'r');
hold on
xlabel('Time');
ylabel('State Vector');
legend('true state 5','Estimated state 5','true state 6','Estimated state 6');
title('State #5 and State #6 Vs Time plot');

figure,
plot(X(7,:),'k');
hold on
plot(X_cap(7,:),'b');
hold on
plot(X(8,:),'g');
hold on
plot(X_cap(8,:),'r');
hold on
xlabel('Time');
ylabel('State Vector');
legend('true state 7','Estimated state 7','true state 8','Estimated state 8');
title('State #7 and State #8 Vs Time plot');

figure,
plot(X(9,:),'k');
hold on
plot(X_cap(9,:),'b');
hold on
plot(X(10,:),'g');
hold on
plot(X_cap(10,:),'r');
hold on
xlabel('Time');
ylabel('State Vector');
legend('true state 9','Estimated state 9','true state 10','Estimated state 10');
title('State #9 and State #10 Vs Time plot');
 

% plot of the ture output , mesured output, and estimated output


figure,
plot(Y(1,:),'r');
hold on
plot(Z(1,:),'b');
hold on
plot(Y_cap(1,:),'g');
hold on
xlabel('Time');
ylabel('Output Vector');
legend('True output','Measured Output','Estimated Output');
title('Output Vs Time plot for State 1');

figure,
plot(Y(2,:),'r');
hold on
plot(Z(2,:),'b');
hold on
plot(Y_cap(2,:),'g');
hold on
xlabel('Time');
ylabel('Output Vector');
legend('True output','Measured Output','Estimated Output');
title('Output Vs Time plot for State 2');

figure,
plot(Y(3,:),'r');
hold on
plot(Z(3,:),'b');
hold on
plot(Y_cap(3,:),'g');
hold on
xlabel('Time');
ylabel('Output Vector');
legend('True output','Measured Output','Estimated Output');
title('Output Vs Time plot for State 3');

 % plot of  norm of error covarince matrix kalman gain matrix
 
stdystate_ndx = round( 0.030/sys.ts );
stdystate_cov_est = var( X_cap(:,stdystate_ndx:end)-X(:,stdystate_ndx:end), 0, 2 );
approx_norm_stdystate_cov = max( stdystate_cov_est );
figure,
semilogy(norm_P,'r');grid
hold on
semilogy(norm_K,'b');grid
hold on
semilogy(stdystate_ndx:20000,approx_norm_stdystate_cov,'g');grid
hold on
xlabel('Time');
ylabel('Error Covariance and Kalman Gain in logscale');
legend('Norm of Error Covariance','Norm of Kalman Gain','State Estimation Error');
title('Error Covariance and Kalman Gain Vs Time plot');





    
