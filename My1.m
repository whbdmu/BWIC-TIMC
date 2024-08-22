function res=My1(S,G,truth,c,ind_folds,p,mode,r)
%% Initialization
[~,n2,~]=size(S);
for i=1:n2
    S_sss(:,:,i)=G{i}'*S{i}*G{i};
end
numInst = length(truth);
dim=size(S_sss);
[n1,n2,n3]=size(S_sss);
mu=1e-4;
tol = 1e-10; max_iter = 200; rho = 1.1;
Q1=zeros(dim);
A=zeros(n2); 
T=zeros(dim); 
iter=0;
alpha=1/n3*ones(1,n3);
Z=S_sss;
beta=ones(n1,1);
alpha_r = alpha.^r;
tempZ = zeros(n1);  
for i=1:max_iter
    iter=iter+1;
    X_k=Z;
    Z_k=T;
    loss=0;
    %% ------Update Z------%%
   
    for iv = 1:length(S)
        W = ones(numInst,numInst);
        ind_0 = find(ind_folds(:,iv) == 0);  % indexes of misssing instances
        W(:,ind_0) = 0;
        W(ind_0,:) = 0;
        for j=1:n3
            Linshi_L = (alpha_r(iv)*S_sss(:,:,j).*W+alpha_r(iv)*A+0.5*mu*(T(:,:,j)-(1/mu)*Q1(:,:,j)))./(alpha_r(iv)+alpha_r(iv)*W+0.5*mu);
        for num = 1:numInst
            indnum = [1:numInst];
            indnum(num) = [];
            Z(indnum',num,j) = (EProjSimplex_new(Linshi_L(indnum',num)'))';
        end
        tempZ = tempZ + alpha_r(iv) * Z(:,:,j);
        end
        clear Linshi_L 
    end

    
    %% ------Update  A------%%
    tempZ = tempZ ./ sum(alpha_r);
    A = tempZ - diag(diag(tempZ));
    A = max(0.5 * (A + A'), 0);
   
    %% ------Update  T------%%
    [T,~,~] = prox_tnn(Z+Q1/mu,beta/mu,p,mode);
   
    %% ------Update F------%%
    temp_W=A;
    temp_W = (temp_W+temp_W')/2;                                                       
    L_D = diag(sum(temp_W));
    L_Z = L_D-temp_W;
    [F, ~, ~]=eig1(L_Z, c, 0);
    
 
    
    %% ------Update alpha------%%
    Rec_error = zeros(1,n3);
    for iv = 1:n3
      Rec_error(iv) = norm((Z(:,:,iv)-A),'fro')^2+norm((Z(:,:,iv)-S_sss(:,:,iv)).*W,'fro')^2;
    end
    H = bsxfun(@power,Rec_error, 1/(1-r));     
    alpha = bsxfun(@rdivide,H,sum(H)); 
    alpha_r = alpha.^r; 
   clear H
    %% ------Checking Convergence-------%%
    chgZ=max(abs(T(:)-Z_k(:)));
    chgX_Z=max(abs(Z(:)-T(:)));
    chg=max([ chgZ chgX_Z]);
    for v = 1:n3
         loss=loss+alpha_r(iv)*norm((Z(:,:,iv)-S_sss(:,:,iv)).*W,'fro')^2+alpha_r(iv)*norm(Z(:,:,iv)-A,'fro')^2;
    end
    loss=loss+trace(F'*L_Z*F);
    if iter == 1 || mod(iter, 10) == 0
        disp(['iter ' num2str(iter) ', mu = ' num2str(mu) ', chg = ' num2str(chg) ',  chgZ = ' num2str(chgZ) ',chgX_Z = ' num2str(chgX_Z) ]);
    end

    %% -------Update Lagrange multiplier-------%%
   
    Q1=Q1+0.5*mu*(Z-T);
    mu=rho*mu;

    %% --------Clustering-------%%
    new_F = F;
    norm_mat = repmat(sqrt(sum(new_F.*new_F,2)),1,size(new_F,2));
    for i = 1:size(norm_mat,1)
        if (norm_mat(i,1)==0)
            norm_mat(i,:) = 1;
        end
    end
    new_F = new_F./norm_mat;
    repeat = 5;
    for iter_c = 1:repeat
        pre_labels    = kmeans(real(new_F),c,'emptyaction','singleton','replicates',20,'display','off');
        result_LatLRR = ClusteringMeasure(truth, pre_labels);
        AC(iter_c)    = result_LatLRR(1)*100;
        MIhat(iter_c) = result_LatLRR(2)*100;
        Purity(iter_c)= result_LatLRR(3)*100;
    end
    mean_ACC = mean(AC);
    mean_NMI = mean(MIhat);
    mean_PUR = mean(Purity);
    res=[mean_ACC mean_NMI mean_PUR];
    if iter == 1 || mean_ACC>40
        disp(['iter ' num2str(iter) ', mean_ACC = ' num2str(mean_ACC) ', mean_NMI = ' num2str(mean_NMI) ', mean_PUR = ' num2str(mean_PUR)]);
        disp(['---------------------------------------------------------------------------------------'])
    end
end