
%random a directed ER network
%rand(n,p)  and rand_ER(n,p,E) ignore p 
function adj = rand_ER(n,p,E)

 

switch nargin
  
    case 2 % the number of nodes and the probability of attachment, n, p
   adj= rand(n,n) < p;
   adj=double(adj);
    for i=1:n
       adj(i,i)=0;
    end 
        

  
    case 3 % fixed number of nodes and edges, n, E   
        adj=zeros(n); % initialize adjacency matrix
        while numedges(adj) < E
            i=randi(n); j=randi(n);
            if i==j | adj(i,j)>0; continue; end  % do not allow self-loops or double edges
            adj(i,j)=1; 
        end
    
    
end  % end nargin 