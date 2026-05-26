function p=partition(n,k)
% p=partition(n,k)
%
% given n and k returns a vector p=[p1, p2, ... p_k] with random integers
% from 1 to (n-(k-1)) such that sum(p)=n
%
% This solution was found by the student Paulo Rodrigues (2223323-not confirmed!) 26May2023
% with some adjustments to make it consistent and completely vectorized
%p=diff([0, sort(randi(n-k+1,1,1,k-1)),n]);
%p=diff([0, sort(randperm(n-k+1,1,k-1)), n]);

p=diff([0 sort(randperm(n-1,k-1)) n]);