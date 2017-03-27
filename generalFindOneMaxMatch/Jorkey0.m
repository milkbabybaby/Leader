clear all;
clc;
load('A.mat')
[V,J]=jordan(A); %V-1AV=J; Q=v-1B
