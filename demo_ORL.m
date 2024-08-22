clear;
clc;
addpath(genpath('Measure/'));
addpath(genpath('tSVD/'));
%Dataname = 'ORL';

%% ------INCOMPLETE DATA PREPARATION------ %%
load('dataset\ORL_4views.mat'); truth=truth';
load('dataset\ORL_fold_example.mat');
ind_folds =A; %index matrix
clear A
numClust = length(unique(truth));
num_view = length(X);
[numFold,numInst]=size(ind_folds);
c=numClust;
res = [];
rand('seed',931316785); %You need to choose the appropriate random seed  
%% --------GRAPH CONSTRUCT-------- %%
for iv = 1:num_view
    X1 = X{iv}';
    X1 = NormalizeFea(X1,1);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = [];
    Y{iv} = X1';
    W1 = eye(size(ind_folds,1));
    W1(ind_0,:) = [];
    G{iv} = W1;
end
S_temp=graph_construction(Y);

%% ----------TRAIN-----------%%
for r =1.5
    for p= 0.4
            mode=2;
            [res] = My1(S_temp,G,truth,c,ind_folds,p,mode,r);
            fprintf(fid,'r: %f ',r);
            fprintf(fid,'p: %f ',p);
           fprintf(fid,'res1: %g %g %g  \n',res(:,1:3));
    end
end
