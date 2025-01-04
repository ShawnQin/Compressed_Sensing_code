function S = geneSimiMatr(N, a, b)
% this function reteturn a N by N similarity matrix, whose elements
% represents simlarity between two odors
% S is symmetrical, and Sij ~ [0,1], generated from a Beta distributon
%
% S     NxN symmetric matrix, with diagnoal as 0
% a     parameter of gamma distribution
% b     parameter of gamma distribution
% X = betarnd(a,b,[N*(N-1)/2,1]);  %upper diagnal elements
S = tril(2*betarnd(a,b,[N,N]) - 1,-1);   %initialize the similarity matrix
S = S'+ S;