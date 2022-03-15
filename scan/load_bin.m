function []=load_bin(channels)
% Load raw data from .bin file

nPacks=100;
rawVectorLength=nPacks*216;
nSamps=nPacks*3;

%%% Spike data %%%
% Create reference matrix to pull spike sample data out of raw file vector %
a=repmat(1:216:nPacks*216,1,3);
a=sort(a);
b=repmat([0 64 128],1,nPacks);
chRef=a+b;

