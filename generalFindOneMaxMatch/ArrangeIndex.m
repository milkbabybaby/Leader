function [ row,column ] = ArrangeIndex( rowFake,columnFake,counter )
%���ھ������޽�����ÿ�εõ������кŶ���Ҫ��������
%   Detailed explanation goes here

%% ��������
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

