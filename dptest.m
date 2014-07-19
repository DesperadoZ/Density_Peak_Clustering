%% density peak test 1 
%----------------------------------------------------------------------
% This script is a simple MATLAB implementation of 
% A.Rodriguez and  A.Laio's 'density peak finding clustering' algorithm[1]
% 
% Class 'nodez' is required to store structured sample/observations
% The input dataset should be a m by n feature matrix which rows are
% observations/samples and columns are features
% The dataset for import must contain only 1 feature matrix to avoid
% confusion. 
%
% [1] Rodriguez, A. and A. Laio (2014). "Clustering by fast search and  
% find of density peaks." Science 344(6191): 1492-1496.
% Last modifed by: zhxngzhx on 18-07-2014
% -------------------zhxngzhx2014--------------------------------------
%% Initialization

clear all; clc; close all

%-------- Change this section for different dataset ----------------------
%load dptest_dataset2

s=uiimport;
ss=struct2cell(s);
FeatMat=ss{1};

options.Resize='on';
defaultanswer={'10','1','1'};
X=inputdlg({'Number of features','Number of labels',...
    'Local distance reference(%)'},'Import parameters',...
    1,defaultanswer);

N_fcol=str2num(X{1});  % Number of features
N_lcol=str2num(X{2});   % Number of labels
pn=str2num(X{3});        % define 'local distance' in % of sample size
FeatMat=FeatMat(:,1:(N_fcol+N_lcol));     % Import dataset with or without label as the last column(s)
[m,n]=size(FeatMat);

%% Input Normalization
normX=questdlg('Would you like to normalized the matrix first?',...
    'Data normalization','Sure,why not','Nah, cant be bothered',...
    'Sure,why not');
switch normX
    case 'Sure,why not'
        FeatMat_N=FeatMat(:,1:N_fcol)-repmat(min(FeatMat(:,1:N_fcol)),[m,1]);
        FeatMat_N=FeatMat(:,1:N_fcol)./repmat((max(FeatMat(:,1:N_fcol))-...
            min(FeatMat(:,1:N_fcol))),[m,1]);
        FeatMat(:,1:N_fcol)=FeatMat_N;
    case 'Nah, cant be bothered'
        warndlg('Normalization may greatly alter the performance','Warn1'...
            ,'modal');
        
end
        



%-------------------------------------------------------------------------

%% Main


%---------------------------Manual Tuning---------------------------------
%pn=1;                       % define 'local distance' in % of sample size
%-------------------------------------------------------------------------

Dmat=zeros(m);              % Distance/weight matrix initialization
Dvec=zeros(sum(1:(m-1)),1); % Distance/weight vector
dpnode=nodez.empty;         % Create empty nodez object array
nclass=0;                   % Number of clusters
cmap=colormap;


j=1;                   % loop index for Dvec

wb=waitbar(0,'Updating distances...');
for i=1:m
    
    dpnode(i).index=i;
    dpnode(i).id='TESTING DATA NO DESCRIPTION';
    dpnode(i).feat=FeatMat(i,1:N_fcol); % select the feature cols
    dpnode(i).label_ref=FeatMat(i,(N_fcol+1):(N_fcol+N_lcol));  % select label col(s), ignore this property if no label column in the dataset
    
    dpnode(i).up_distz(FeatMat(:,1:N_fcol)); % update distances   
    
    Dmat(i,:)=dpnode(i).distz;  % construct distance matrix
      
    if i<m
        tempv=dpnode(i).distz(i+1:end);
        nt=length(tempv);
        Dvec(j:(j+(nt-1)))=tempv; % construct distance vector (non-zero elements)
        j=j+nt;
    end
    
    waitbar(i/m,wb);
end 
close(wb);

dc=prctile(sort(Dvec),pn);  % calculate cutoff distance

wb=waitbar(0,'Updating rho...');
for i=1:m
    dpnode(i).up_rho(dc);   % updates local density
    waitbar(i/m,wb);
end
close(wb);
wb=waitbar(0,'Updating delta...');
for i=1:m
    dpnode(i).up_delta(dpnode); % updates delta distance
    if isinf(dpnode(i).delta);
        dpnode(i).delta=max(Dvec);  % special case when the node has the largest rho
    end
    waitbar(i/m,wb);
end
close(wb);

rho=[dpnode.rho];
delta=[dpnode.delta];

figure(1)
subplot(2,1,1)
p1=plot(rho(:),delta(:),'o','Markersize',5,'MarkerFaceColor','k');
title('Decision Graph: Please select cluster centres','FontSize',15.0);
xlabel('\rho Local density');
ylabel('\delta Distance Delta');
hold on;


rect=getrect(1); % request user input for finding cluster centres (lower bounds for rho and delta)
rhomin=rect(1);
deltamin=rect(2);


wb=waitbar(0,'Finding cluster centres...');
for i=1:m
    if (dpnode(i).rho>rhomin) && (dpnode(i).delta>deltamin)
        nclass=nclass+1;
        dpnode(i).label=nclass; % label cluster centres
    end
    waitbar(i/m,wb);
end
close(wb);

[rho_d,ordrho]=sort(rho,'descend');

wb=waitbar(0,'Clustering...');
for i=1:m
    ii=ordrho(i);
    if isempty(dpnode(ii).label)
        dpnode(ii).label=dpnode(dpnode(ii).d_neigh).label; % label the rest of the nodes based on the class of the nearest higher density node
    end
    waitbar(i/m,wb);
end
close(wb);


class=[dpnode.label];
D_out=zeros(nclass,1);

wb=waitbar(0,'Calculating density threshold...');
for i=1:m                       % find density threshold
    for ii=1:m
        if (class(i)~=class(ii)) && Dmat(i,ii)<=dc
            avg_rho=(rho(i)+rho(ii))/2;
            if avg_rho>D_out(class(i))
                D_out(class(i))=avg_rho;
            end
            if avg_rho>D_out(class(ii))
                D_out(class(ii))=avg_rho;
            end
        end
    end
    waitbar(i/m,wb);
end
close(wb);

wb=waitbar(0,'Locating outliers...');

for i=i:m
    if(rho(i)<D_out(class(i)))
        dpnode(i).label=0;      % 'halo' identified due to low density
    end
    waitbar(i/m,wb);
end
close(wb);

for i=1:nclass 
    hh=findobj(dpnode,'label',i);
    hindex=[hh.index];
    assignin('base',['class',num2str(i)],hindex); % save index vector for each class
    
end




for i=1:nclass  % Redraw fig.1a to illustrate clustering
    c=round((i*64)/nclass);
    subplot(2,1,1);
    P=eval(['class',num2str(i)]);
    
    plot(rho(P),delta(P),'o','MarkerSize',8,'MarkerFaceColor',...
        cmap(c,:),'MarkerEdgeColor',cmap(c-5,:));
end

wb=waitbar(0,'Visualizing output please wait...');
subplot(2,1,2)
Y1=mdscale(Dmat,2,'criterion','metricstress');
plot(Y1(:,1),Y1(:,2),'*','MarkerSize',5,'MarkerFaceColor','k',...
    'MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')

for i=1:nclass % Redraw fig.1b to illustrate clustering
    c=round((i*64)/nclass);
    subplot(2,1,2);
    hold on;
    P=eval(['class',num2str(i)]);
    
    plot(Y1(P,1),Y1(P,2),'o','MarkerSize',8,'MarkerFaceColor',...
        cmap(c,:),'MarkerEdgeColor',cmap(c-5,:));
    
end
% gname;
waitbar(1,wb);
close(wb);

msg=msgbox([num2str(nclass),' clusters have been found!'],'Summary');


            
            
            






    








    


