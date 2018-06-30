function [max_corr_valid,corr_CAI_valid,best,jump,corr_CAI_test,corr_test] = greedy_cor_better(A,PA1)
%% 
L=length(A); %number of genes
%creating random sets
rando = randperm(L);
train = rando(1,1: ceil(L/3));
test = rando(1,ceil(L/3)+1 : floor(2*L/3));
valid = rando(1,ceil(2*L/3): L) ;
%
CAI = A(:,1); % CAI should be the first column in the array A
best = []; % initializing traceback vector
jump = []; % initializing 'difference in correlation due to extra parameters' vector (funny name)
% finding the number of colums
b=size(A); 
colums = b(2); 
%finding length of the randomized sets
L1 = length(PA1(train));
L2 = length(PA1(test)) ;
L3 = length(PA1(valid));
%%
mx_vec = [];
for i = 2:colums
    % finding the current best correlation for comparison with next
    % parameters
    x_train=[ones(L1,1) A(train,best) CAI(train)];
    [beta]=regress(x_train,PA1(train));
    x_test = [ones(L2,1) A(test,best) CAI(test)];
    Y= x_test*beta;
    [max_corr,max_Pval]=corr(transpose(PA1(test))',transpose(Y)','Type','Spearman','Rows','complete');
    max_corr_last = max_corr ;
    % 
    max_index = [] ; % initializing a new max parameter's index for every iteration
    mx_vec= [mx_vec,max_corr];
    for j = 2:colums
        if any(best == j) == 0
            % computing the temporary correlation with the j-th parameter
            x_temp_train=[ones(L1,1) A(train,best) CAI(train) A(train,j)];
            x_temp_test=[ones(L2,1) A(test,best) CAI(test) A(test,j)];
            [row,col] = size(x_temp_train);
            if rank(x_temp_train) < col || rank(x_temp_test) < col 
                continue
            else
                
                [row,col]=size(x_temp_train);
                [beta_temp]=regress(x_temp_train,PA1(train));
                Y_temp = x_temp_test*beta_temp;
                [temp_corrPA,temp_PvalPA]=corr(transpose(PA1(test))',transpose(Y_temp)','Type','Spearman','Rows','complete');
            end
           
            %
            if temp_corrPA > max_corr
                max_corr = temp_corrPA;
                max_index = j ;
            end
        end
    end
    best = [best,max_index];
    jump = [jump,max_corr - max_corr_last];
end

% computing the correlation with all the parameters that we found on the
% validation group and compare it to only using the CAI parameter
x_max_train=[ones(L1,1) A(train,best) CAI(train)];
[beta_best]=regress(x_max_train,PA1(train);
x_max_valid=[ones(L3,1) A(valid,best) CAI(valid)];
Y_best = x_max_valid*beta_best;
[max_corr_valid,max_Pval]=corr(transpose(PA1(valid))',transpose(Y_best)','Type','Spearman','Rows','complete');

[corr_CAI_valid,PVAL_CAI]=corr(transpose(PA1(valid))',transpose(CAI(valid))','Type','Spearman','Rows','complete');
%
corr_CAI_test=mx_vec(1);
corr_test=mx_vec(end);

figure; 
scatter(PA1(valid),Y_best,'r')
title('Predicted Values as a function of Protein Abundance ')
set(gca,'xscale','log')
xlabel('log(PA)')
ylabel('Predictor')
figure; 
scatter(PA1(valid),CAI(valid),'b')
title('CAI Values as a function of Protein Abundance ')
set(gca,'xscale','log')
xlabel('log(PA)')
ylabel('CAI')

end

