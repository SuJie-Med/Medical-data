clc;
clear;
close all;
%% set the parameters
K=4;%%% the number of elements in a cluster
a=1;b=0;
%read the data
A0 = fopen('data60.csv');
data60 = textscan(A0, ' %d32 %s %d32 %d32 %s %s','delimiter', ',')
fclose(A0);
%% initial
X.N=data60{1,3};
X.X0=data60{1,2};
X.X1=data60{1,4};
X.X2=data60{1,5};
X.X3=data60{1,6};
%X.X4=data60{1,7};
nn=size(X.N,1);%%%nodes in graph
%% read the network G0
A=ones(nn,nn)*1000;
% A=zeros(nn,nn)*1000;
%%%draw graph 
%%%%%%
%% draw graph
% G=graph(A)
% plot(G);
%%
G0=A(1:nn,1:nn);
X0.N=X.N(1:nn);
X0.X0=X.X0(1:nn);
X0.X1=X.X1(1:nn);
X0.X2=X.X2(1:nn);
X0.X3=X.X3(1:nn);
% X0.X4=X.X4(1:nn);
G1=A;
%% define attributes
attr.c0={'F','M'};
attr.c1={'32***','33***','34***','72***','75***'};
attr.c2={'M','S','D'};
attr.c3={'Y','N'};
%attr.c4={'Y','N'};
%attr.N=X.N;
attr.N=X0.N;
%[GilS,NSilS,S0]=ComputeLoss(G0,attr,X0,K,a,b);
NN=nn/2;

%%
for K=2:2:floor(NN)  
  for a=1:2:11
    [GilS(K,a),NSilS(K,a),S0]=ComputeLoss(G1,attr,X0,K,(a-1)/10,1-(a-1)/10);
     SCluster{K,a}=S0;
  end; 
end;

 save('SCluster.mat','SCluster');
 save('GilS.mat','GilS');
 save('NSilS.mat','NSilS');


