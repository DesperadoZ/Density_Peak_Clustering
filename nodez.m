classdef nodez<handle
    %----------------------------------------------------------------------
    % Nodez objects can be used to repersent observations from any dataset
    % This particular version was modified to suit A.Rodriguez and
    % A.Laio's 'density peak finding clustering' algorithm [1]
    %
    % Please note that if the node has highest density, the delta value
    % will be marked as 'inf'
    % 
    %
    % [1] Rodriguez, A. and A. Laio (2014). "Clustering by fast search and  
    % find of density peaks." Science 344(6191): 1492-1496.
    % Last modifed by: ZhxngZhx on 18-07-2014 
    % -------------------ZhxngZhx2014--------------------------------------
    
    properties
        index=0;    
        id='';
        feat=[];    %feature vector
        label_ref=[];   % training/actual label 
        label=[];   %clustering label
        distz=[];   %distance/weight vector 
        rho=0;      %local density
        delta=0;    %distance to nearest higher density node
        
        
       
        
    end
    
    properties (Dependent = true, SetAccess = private)
        dim
        d_neigh=0;  %the index of the nearest higher density node
        
    end
    
    methods
        function dim=get.dim(obj)
            dim=length(obj.feat);
        end
        
        function set.feat(obj,val)
            obj.feat=val(:)';
        end
        
        function set.label(obj,val)
            obj.label=val(:)';
        end
        
        function up_distz(obj,feat_matrix)
            [m,n]=size(feat_matrix);
            if obj.dim~=n
                error('The number of features must agree!')
            end
            
            for i=1:m
                obj.distz(i)=sqrt(sum((obj.feat-feat_matrix(i,:)).^2));
            end
        end
        
        function up_rho(obj,dc)
            n=obj.dim;
            obj.rho=0;
            for i=1:n
                obj.rho=obj.rho+exp(-(obj.distz(i)/dc)^2);
            end
        end
        
        function up_delta(obj,objarr)
            rhov=[objarr.rho];
            pos=find(rhov>obj.rho);
            if isempty(pos)
                obj.delta=inf;
            else
            obj.delta=min(obj.distz(pos));
            end
        end
        
        function d_neigh=get.d_neigh(obj)
            
            if ~isinf(obj.delta)
                  d_neigh=find(obj.distz==obj.delta);
            else
                  d_neigh=obj.index;
            end
        end
            
            
            
        

    end
    
end

