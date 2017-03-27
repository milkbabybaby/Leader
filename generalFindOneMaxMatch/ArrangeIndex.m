function [ row,column ] = ArrangeIndex( rowFake,columnFake,counter )
%由于矩阵在修建，故每次得到的行列号都需要进行修正
%   Detailed explanation goes here

%% 函数主体
for i = counter : -1 : 2
    for j = 1 : 1 : i-1
        if rowFake(i) >= rowFake(i-j)
            rowFake(i) = rowFake(i) + 1;
        end
        if columnFake(i) >= columnFake(i-j)
            columnFake(i) = columnFake(i) + 1;
        end
    end
end

row = rowFake(1:counter);
column = columnFake(1:counter);

end

