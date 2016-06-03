A= [-1 250 0 0 0 0 0;...
    -250 -1 0 0 0 0 0;...
    0 0 -1 190 0 0 0;...
    0 0 -190 -1 0 0 0;...
    0 0 0 0 -1 150 0;...
    0 0 0 0 -150 -1 0;...
    0 0 0 0 0 0 -300];
B= [21.1 10.2;...
    7.63 12.7;...
    7.65 4.21;...
    2.17 8.32;...
    6.45 0.85;...
    1.52 5.13;...
    25 21.3]; 
C= [7.63 21.1 2.17 7.65 1.52 6.45 25];
D= [0 0];
n=size(A); 
olpoles=eig(A);
clpoles=[-50-250i...
         -50+250i...
         -38-190i...
         -38+190i...
         -30-150i...
         -30+150i...
            -300];
disp('the open loop pole loaction are:');
disp(olpoles);
disp('the desired closed loop pole loaction are:');
disp(clpoles);
olcharcoff=poly(olpoles);
clcharcoff=poly(clpoles);
%loop for singular value decomposition

psi=zeros(9,7);
 for i=1:n
    p=( cat(2,( clpoles(i)*eye(n)-A ),B));
    [U,S,V]=svd(p);
    psi(:,i)=V(:,8)+V(:,9);

 end

 disp('psi=')
disp(psi)
disp('the gain matrix K=');
K=psi(8:end,:)/( psi(1:7,:) );

norm_K=norm(K);
disp('the norm of the gain matrix using mimo null space method is');
disp(norm_K);
disp('the norm of this gain matrix is comaparitively less than what achived ');
disp('by the method of ackman"s formula hence the strenght of the signal to');
disp('achive the gain in this method is less ,hence benificial');


%verification of the gain matrix

cleig=eig(A-B*K);
disp('the closed loop eigon values are');
disp(cleig);


 

 