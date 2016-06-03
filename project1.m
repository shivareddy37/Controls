%program to compute the state trastion matrix
prompt=('enter matrix A of the state equation A=  ');
A=input(prompt);
disp(A);
prompt=('enter matrix B of the state equation B=  ');
B=input(prompt);
disp(B);
prompt=('enter matrix C of the state equation C=  ');
C=input(prompt);
disp(C);
prompt=('enter matrix D of the state equation D=  ');
D=input(prompt);
disp(D);
sys=ss(A,B,C,D);
%disp(sys);
n=size(A);
Qc=zeros(n);
Qc(:,1)=B;
for k=2:n
    Qc(:,k)=A*Qc(:,k-1);
end
p=poly(A);
disp('coffcients of characterstic polynomial are p =');disp(p);
r=roots(p);
disp('roots of the characterstic polynomial are r =');disp(r);
c=[1 0 0 0 0 0 0];
N=toeplitz(c,p);
disp(N);
T=N*Qc;
disp('the transformation matrix is T=');disp(T);
disp
