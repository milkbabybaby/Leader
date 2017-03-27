function [b]=scalefree(ni,m,de)
b=zeros(m,1);
dp=zeros(ni+1,1);
sk=0;
for i=1:ni
    sk=sk+de(i); %sk 总度数
end
dp(1,1)=0;
for i=1:ni
    dp(i+1,1)=dp(i,1)+de(i,1);
end %dp=[0,dk]'
is=1;
while is<=m 
r=rand(1,1);
r=fix(r*sk+1);  %在[1,sk]按平均分布 随机选择一个数
    for i=1:ni
        if r>dp(i,1)&r<=dp(i+1,1)
            it=i;
        end
    end
    pd=0;  
    for j=1:is
        if it==b(j,1)
            pd=1;
        end
    end  %查重;d=1 it已经在b里了，在来一遍 否则就要赋值
    if pd==0
       b(is,1)=it;
       is=is+1;
    else
        is=is;
    end
end

