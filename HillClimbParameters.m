clear; clc;
% %% importing the data
[~,~,PT]=xlsread('C:\Users\mishael\Desktop\matlab files\Project\Ptable'); 
[np,~]=size(PT); %number of proteins
load('parameters_new.mat');
% 
% %% double AA table
% [~,~,DT]=xlsread('C:\Users\mishael\Desktop\studys\project\II\ST','sheet1'); 
% [nd,~]=size(DT); %size of excel sheet.
% ss2p=cell2mat(DT(2:nd,2));
% ss2n=find(ss2p<=0.001);
% ss2n=ss2n+1;
% ss2c=DT(ss2n,1);
% ss2n=ss2n-1;
% ss2p=ss2p(ss2n);
% ss2c{length(ss2c),3}=[]; % adding more vector to be used later
% ss2c(:,2:3)={zeros(np,1)};
% 
% %% triple AA table
% [~,~,TT]=xlsread('C:\Users\mishael\Desktop\studys\project\II\ST','sheet2'); 
% [nt,~]=size(TT); %size of excel sheet.
% ss3p=cell2mat(TT(2:nt,2));
% ss3n=find(ss3p<=0.001);
% ss3n=ss3n+1;
% ss3c=TT(ss3n,1);
% ss3n=ss3n-1;
% ss3p=ss3p(ss3n);
% ss3c{length(ss3c),3}=[];
% ss3c(:,2:3)={zeros(np,1)};
% %% finding the parameters
% for i=2:1:np
%     P=nt2aa(PT{i,7});% the translated protein for each gene
%     if length(P)>49
%         for j =1:1:49 %first 50 dipeptides count
%             AA=P(j:j+1);
%             IndexC = strfind(ss2c(:,1), AA);
%             Pnum = find(not(cellfun('isempty', IndexC)));
%                 if isempty(Pnum)==0
%                     ss2c{Pnum,2}(i)=ss2c{Pnum,2}(i)+1;
%                 end
%         end
%     end
%     for j =50:1:length(P)-1 %after 50 dipeptides
%         AA=P(j:j+1);
%         IndexC = strfind(ss2c(:,1), AA);
%         Pnum = find(not(cellfun('isempty', IndexC)));
%             if isempty(Pnum)==0
%                 ss2c{Pnum,3}(i)=ss2c{Pnum,3}(i)+1;
%             end
%     end
%     if length(P)>49
%         for j =1:1:48 % first 50 tripeptides count
%             AA=P(j:j+2);
%             IndexC = strfind(ss3c(:,1), AA);
%             Pnum = find(not(cellfun('isempty', IndexC)));
%                 if isempty(Pnum)==0
%                     ss3c{Pnum,2}(i)=ss3c{Pnum,2}(i)+1;
%                 end
%         end  
%     end
%     for j =49:1:length(P)-2 % last peptides count
%         AA=P(j:j+2);
%         IndexC = strfind(ss3c(:,1), AA);
%         Pnum = find(not(cellfun('isempty', IndexC)));
%             if isempty(Pnum)==0
%                 ss3c{Pnum,3}(i)=ss3c{Pnum,3}(i)+1;
%             end
%     end
%     for j=1:1:53      
%         ss2c{j,3}(i)=ss2c{j,3}(i)/length(P);
%     end
%     for j=1:1:247
%         ss3c{j,3}(i)=ss3c{j,3}(i)/length(P);
%     end
% end
% 
% % making a cell array of all the data
% 
% %adding the number of appearence
% m=53*2+247*2;% size of new cell array
% parameters=cell(np,m); % the parameters cell array
% counter=0;% counts where we stoped the loop
% for i=1:1:53
%     parameters{1,i}= strcat(ss2c{i-counter,1},'_first');
%     for j=2:1:np
%         parameters{j,i}= ss2c{i-counter,2}(j-1);
%     end
% end
% counter=i;
% for i=54:1:53+counter
%     parameters{1,i}= strcat(ss2c{i-counter,1},'_last');
%     for j=2:1:np
%         parameters{j,i}= ss2c{i-counter,3}(j-1);
%     end
% end
% counter=i;
% for i=1+counter:1:247+counter
%     parameters{1,i}= strcat(ss3c{i-counter,1},'_first');
%     for j=2:1:np
%         parameters{j,i}= ss3c{i-counter,2}(j-1);
%     end
% end
% counter=i;
% for i=1+counter:1:247+counter
%     parameters{1,i}= strcat(ss3c{i-counter,1},'3_last');
%     for j=2:1:np
%         parameters{j,i}= ss3c{i-counter,3}(j-1);
%     end
% end
% %GC
% [n,m]=size(parameters);
% parameters{1,m+1}='GC';
% for i =2:1:n
%     parameters{i,m+1}=GC(i-1);
% end
% %Len
% [n,m]=size(parameters);
% parameters{1,m+1}='length';
% for i =2:1:n
%     parameters{i,m+1}=LEN(i-1);
% end
% %A_p
% [n,m]=size(parameters);
% parameters{1,m+1}='A_p';
% for i =2:1:n
%     parameters{i,m+1}=A_P(i-1);
% end
% %AT
% [n,m]=size(parameters);
% parameters{1,m+1}='AT';
% for i =2:1:n
%     parameters{i,m+1}=AT(i-1);
% end
% %C_P
% [n,m]=size(parameters);
% parameters{1,m+1}='C_P';
% for i =2:1:n
%     parameters{i,m+1}=C_P(i-1);
% end
% %G_P
% [n,m]=size(parameters);
% parameters{1,m+1}='G_P';
% for i =2:1:n
%     parameters{i,m+1}=G_P(i-1);
% end
% %the first 50 II score and last 50, normal way
% [n,m]=size(parameters);
% parameters{1,m+1}='II2_first50_normal';
% for i =2:1:n
%     parameters{i,m+1}=II21(i-1,1);
% end
% [n,m]=size(parameters);
% parameters{1,m+1}='II2_after50_normal';
% for i =2:1:n
%     parameters{i,m+1}=II21(i-1,2);
% end
% 
% %same but now with median
% [n,m]=size(parameters);
% parameters{1,m+1}='II2_first50_med';
% for i =2:1:n
%     parameters{i,m+1}=mat_median_II2(i-1,1);
% end
% [n,m]=size(parameters);
% parameters{1,m+1}='II2_after50_med';
% for i =2:1:n
%     parameters{i,m+1}=mat_median_II2(i-1,2);
% end
% %II2 first #
% [n,m]=size(parameters);
% for i=1:1:10
%     parameters{1,m+i}=strcat('II2_first_',num2str(i));
%     for j =2:1:n
%       parameters{j,m+i}=matfirst(j-1,i);
%     end
% end
% % II2 last #
% [n,m]=size(parameters);
% for i=1:1:50
%     parameters{1,m+i}=strcat('II2_last_',num2str(i));
%     for j =2:1:n
%       parameters{j,m+i}=matlast(j-1,i);
%     end
% end
% % II first first and last normal way
% [n,m]=size(parameters);
% parameters{1,m+1}='II3_first50_normal';
% for i =2:1:n
%     parameters{i,m+1}=II3_mat_normal(i-1,1);
% end
% [n,m]=size(parameters);
% parameters{1,m+1}='II3_after50_normal';
% for i =2:1:n
%     parameters{i,m+1}=II3_mat_normal(i-1,2);
% end
% 
% %now with median
% [n,m]=size(parameters);
% parameters{1,m+1}='II3_first50_med';
% for i =2:1:n
%     parameters{i,m+1}=mat_median_II3(i-1,1);
% end
% [n,m]=size(parameters);
% parameters{1,m+1}='II3_after50_med';
% for i =2:1:n
%     parameters{i,m+1}=mat_median_II3(i-1,2);
% end
% 
% %II3 first #
% [n,m]=size(parameters);
% for i=1:1:10
%     parameters{1,m+i}=strcat('II3_first_',num2str(i));
%     for j =2:1:n
%       parameters{j,m+i}=matfirst_tri(j-1,i);
%     end
% end
% % II3 last #
% [n,m]=size(parameters);
% for i=1:1:50
%     parameters{1,m+i}=strcat('II3_last_',num2str(i));
%     for j =2:1:n
%       parameters{j,m+i}=matlast_tri(j-1,i);
%     end
% end
%CAI
% [n,m]=size(parameters);
% parameters(1:n,m+1)=PT(:,8);
%CAI yoav way

CAI=PT(:,8);
parameters=[CAI parameters_new];
%final results
PA1=PT(2:end,3);
PA2=PT(2:end,4);
PA3=PT(2:end,5);
for i=1:1:length(PA1)
    if strcmp(PA1{i},'nan')
        PA1{i}= NaN;
    end
    if strcmp(PA2{i},'nan')
        PA2{i}= NaN;
    end
    if strcmp(PA3{i},'nan')
        PA3{i}= NaN;
    end
    
end
%removing nans
%PA1
PA1=cell2mat(PA1);
PA2=cell2mat(PA2);
PA3=cell2mat(PA3);
nans=isnan(PA1)+isnan(PA2)+isnan(PA3);
nans(nans>1)=1;
nans=logical(nans);
PA1(nans)=[];
PA2(nans)=[];
PA3(nans)=[];
%dividing by mean
PA1=PA1./mean(PA1);
PA2=PA2./mean(PA2);
PA3=PA3./mean(PA3);
%dividing into 3 subgroups
pameans=(PA1+PA2+PA3)/3;
a=cell2mat(parameters(2:end,:));
a(nans,:)=[];
[max_corr_valid,corr_CAI_valid,best,jump,corr_CAI_test,corr_test] = testcorr(a,pameans);
best_names=cell(1,length(best));
for i =1:1:length(best)
    best_names{i}=parameters{1,best(i)};
end

