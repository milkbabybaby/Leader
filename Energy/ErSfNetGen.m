function [ Net ] = ErSfNetGen( kkk,NumNode,meanDegree )
%产生ER网络和SF(BA）网络
%kkk=1为ER网络，否则为SF（BA）网络
%NumNode为节点数，meanDegree为平均度，即Uin或者Uout
%Net为NumNode*NumNode邻接矩阵

    

    N = NumNode;
    if kkk==1
        P=poissrnd(meanDegree,1,N);
        A=rand(N,N);
        for i=1:N
            p=P(i)/N;
            if p<0
                p=0;
            end
            

            for j=1:N

                
                

                if A(i,j)<p %&& b2>0.5%abs(normrnd(2.5/N,0.001))
                    A(i,j)=1;
                    %             elseif A(i,j)<p  && b2<=0.5
                    %                 A(j,i)=1;
                else
                    A(i,j)=0;
                end
            end

        end
        for i=1:N
            A(i,i)=0;
        end

    else
        n=N;
        a=zeros(n,n);
        m=meanDegree; % The mean degree
        % The initial random network
        n0=21;
        p0=0.25;
        for i=1:n0
            for j=1:n0
                if rand(1,1)<p0
                    a(i,j)=1;
                    %   a(j,i)=1;
                end
            end
        end
        for i=1:n0
            deg(i,1)=sum(a(i,:))+sum(a(:,i));
        end
        %%%%%%%%%%%%%%%%%%%%
        for i=n0:n-1
            b=zeros(m,1);
            [b]=scalefree(i,m,deg);
            for j=1:m
                p=rand;
                if p>=0.5
                    % a(b(j,1),i+1)=1;
                    a(i+1,b(j,1))=1;
                else
                    a(b(j,1),i+1)=1;
                end

                deg(b(j,1),1)=deg(b(j,1),1)+1;
            end
            deg(i+1,1)=m;
        end

        A=a;
        for i=1:N
            A(i,i)=0;
        end

    end
    Net = A;

end

