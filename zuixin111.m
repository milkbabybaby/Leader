clear all
 clc
format long 
%N=100
 
   M=3; N=6
   
A=rand(N,N)
 
% for k=1:N
%     A(k,k)=0
% end
A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);


B=rand(N,M); 
B=I(:,1:M)
% 
% B=[1 0 0 0 0 0; 0 0 0 0 0 1]';
%B=[0 1]';

 %B=[1 0; 0 0; 0 1; 0 0];
I1=eye(M);
tf=1



% WB0=zeros(N,N);
% ot=0.025;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
% end
%  C0=pinv(WB0); 
 
pf=expm(A*tf);

%xf=pf*x0;


% syms s 
% AA=s*I-A;
% invAA=inv(A);
% p=ilaplace(invAA,t);
% p=simplify(p);



iterations=20
v=0.005;



Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  



 
for ii=1:iterations
   
  
    pf=expm(A*tf);
    
 WB0=zeros(N,N);
ot=0.01;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=pinv(WB0); 

   
Cost1(ii)=trace(C*pf*pf');

   
V=eig(A);

 II=eye(N);

for  iii=1:N-1
    
    if V(iii)==conj(V(iii+1)) & V(iii)~=V(iii+1)
        
        II(iii,iii)=0.5;
        
        II(iii,iii+1)=0.5;
           II(iii+1,iii)=0.5*i;
        
        II(iii+1,iii+1)=-0.5*i;
        
  
        end 
end

BB=fliplr(vander(V));




% if rank(diag(V)) < N;
% 
%  break
% else  
BBB=II*BB;

    
CC2=exp(V*tf);

CCC2=II*CC2;
%alpha_tf=pinv(BBB)*CCC2;
alpha_tf=BBB\CCC2;

%alpha_tf =TVreg(BBB,CCC2,100)


Cpf=C*pf;

F111=zeros(N,N);

ot=0.01;
for kk=1:tf/ot
    
CC1=exp(V*ot*kk);
CCC1=II*CC1;
%alpha_t=pinv(BBB)*CCC1;
alpha_t=BBB\CCC1;%




CN1=zeros(N,N);
for kg=1:tf/ot
CN1=expm(A*(tf-kg*ot))*B*B'*expm(A'*(tf-kg*ot))*ot+CN1;
end
CB1=expm(A*tf)-CN1*C*expm(A*tf);







F11=zeros(N,N);
for iiii =1:N-1
    
F1=zeros(N,N);
    for k=1:iiii
        
        F1=A'^(k-1)*CB1*A'^(iiii-k)+F1;

    end

 F11=2*alpha_t (iiii+1)*F1+F11;
end

F111=F11*ot+F111;
  
end 

    
    
  
 
  F2222=zeros(N,N);
for kg=1:tf/ot
     F222=zeros(N,N);
        %CB3=zeros(N,N);
        for kkk=ot:ot:kg*ot
          CB3=CB1*expm(A'*tf)*C*expm(A*(tf-kkk))*B*B';  % CB3=CB1*expm(A'*tf)*C*expm(A*(tf-kkk*ot))*B*B'+CB3;
       CC3=exp(V*(kg*ot-kkk));
CCC3=II*CC3;
alpha_tjiantao=BBB\CCC3;
F22=zeros(N,N);
for iiii=1:N-1
  F2=zeros(N,N);
    for k=1:iiii
        
        F2=A'^(k-1)*CB3*A'^(iiii-k)+F2;
    end

F22=-2*alpha_tjiantao(iiii+1)*F2+F22;
end

F222=F22*ot+F222;
        end


F2222=F222*ot+F2222;
end
 
 

CCG=zeros(N,N);
        for  kk=1:tf/ot  
            CCG=expm(A*(tf-kk*ot))*B*B'*expm(A'*(tf-kk*ot))*ot+CCG;%Ã»¿¼ÂÇ³õÖµÎª0
        end

%CCG=CN1;
   CC4=expm(A'*tf)-expm(A'*tf)*C*CCG;%shangshi        
       % CD3=zeros(N,N);
          F3333=zeros(N,N);
        for kg=1:tf/ot
              F333=zeros(N,N);
        for kkk=ot:ot:kg*ot
            CD3=C*expm(A*tf)*CC4*expm(A*(kg*ot-kkk))*B*B';
        CC5=exp(V*(tf-kkk));
CCC5=II*CC5;
alpha_tfjiantao=BBB\CCC5;
 F33=zeros(N,N);
for iiii=1:N-1
  F3=zeros(N,N);
    for k=1:iiii
        
        F3=A'^(k-1)*CD3*A'^(iiii-k)+F3;
    end

F33=-2*alpha_tfjiantao(iiii+1)*F3+F33;
end


F333=F33*ot+F333;
        end

F3333=F333*ot+F333;
 end
        
 

  
   CE33=zeros(N,N);
        for kk=1:tf/ot
             CE3=zeros(N,N);
        for kkk=ot:ot:kk*ot
            CE3=[expm(A*(tf-kkk))*B*B'*expm(A'*(kk*ot-kkk))*CB1*expm(A'*tf)]*ot+CE3;
        end
        CE33=CE3*ot+CE33;
        end
        F444=zeros(N,N);
for g=1:tf/ot
     F44=zeros(N,N);
for iiii=1:N-1
  F4=zeros(N,N);
    for k=1:iiii
        F4=A'^(k-1)*C*CE33*C*expm(A*g*ot)*B*B'*A'^(iiii-k)+F4;
    end

F44=2*alpha_t(iiii+1)*F4+F44;
end
F444=F44*ot+F444;
end



   
   CF33=zeros(N,N);
        for kk=1:tf/ot
            CF3=zeros(N,N);
        for kkk=ot:ot:kk*ot
            CF3=[expm(A*tf)*CC4*expm(A*(kk*ot-kkk))*B*B'*expm(A'*(tf-kkk))]*ot+CF3;
        end
        CF33=CF3*ot+CF33;
        end
       F555=zeros(N,N);
for g=1:tf/ot
    F55=zeros(N,N);
for iiii=1:N-1
  F5=zeros(N,N);
    for k=1:iiii
        F5=A'^(k-1)*C*CF33*C*expm(A*g*ot)*B*B'*A'^(iiii-k)+F5;
    end

F55=2*alpha_t(iiii+1)*F5+F55;
end
F555=F55*ot+F555;
end
 

 
   F6666=zeros(N,N);
        for kk=1:tf/ot
             F666=zeros(N,N);
     
        for kkk=ot:ot:kk*ot
            CB66=C*expm(A*(tf-kkk))*B*B'*expm(A'*(kk*ot-kkk))*CB1;
      F66=zeros(N,N); 
for iiii=1:N-1
  F6=zeros(N,N);
    for k=1:iiii
        
        F6=A'^(k-1)*CB66*A'^(iiii-k)+F6;
    end

F66=-2*alpha_tf(iiii+1)*F6+F66;
end
F666=F66*ot+F666;
end
F6666=F666*ot+F6666;
 end

        
 F=F3333+F444+F555+F6666+F111+F2222; %F=F111+F2222+F3333+F444+F555+F6666
      
% F=F/norm(F);
% 
% F3=(trace(A'*A)-N)*4*A; % gradient of norm function
% 
% F3=F3/norm(F3);
% 
% F0=reshape((eye(N*N)-F3(:)*pinv(F3(:)))*F(:),N,N); 


%   F=F'+F-diag(diag(F));
%  for kkk=1:N
%     F(kkk,kkk)=0;
% end

A=A-v*F/norm(F);
A=(sqrt((N)/trace(A'*A)))*A; 


ii  
end

% for i=1:length(A)
%     for j=1:length(A)
%        if abs(A(i,j))<1.5*sum(sum(A))/(N*N)
%        A(i,j)=0;
%        end
%     
%     end
% end

figure(1)
  plot((Cost1),'r-*')
  
  figure(2)
  plot(log10(Cost1),'r-*')
  
  
  AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(A);
AAAA-AL

%   id = arrayfun(@(i)[i find(V == conj(V(i)))],[1:numel(V)]','UniformOutput',0);
%  id = cell2mat(id)