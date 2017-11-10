%%Setting up the parameters
alpha=1/3; beta=0.99; sigma=2; 
delta=0.025 ; a_bar=0 ;rho=0.5;
sigma_e= 0.2; l_bar=1;
sigma_z=sigma_e/sqrt(1-rho^2) ; %AR(1) Variance formula
%%
%%Discretize z using Tauchen method
[logz,logzprob]=TAUCHEN(5,rho,sigma_e,3);
z=exp(logz');
Pinv=logzprob^1000;
Ns=sum(Pinv(1,:)*z');

amin=a_bar;
amax=5;
num_a=500;

a=linspace(amin,amax,num_a);

kmin=amin;
kmax=amax;
%%
aggsav=1;
while abs(aggsav)>1e-02;
    k_guess=(kmin+kmax)/2;
    r=alpha*(Ns/k_guess)^(1-alpha)+(1-delta);
    w=(1-alpha)*(k_guess/Ns)^alpha;
    cons=bsxfun(@minus,r*a',a);
    cons=bsxfun(@plus,cons,permute(z*w*l_bar,[1 3 2]));
    ret=cons.^(1-sigma)./(1-sigma);
    ret(cons<0)=-inf;
    
    %Guess value function
    v_guess=zeros(5,num_a);
    
    v_tol=1;
    while v_tol>1e-06;
        v_mat=ret+beta*permute(logzprob*v_guess,[3 2 1]);
        [vfn,pol_idx]=max(v_mat,[],2);
        vfn=permute(vfn,[3 1 2]);
        v_tol=max(max(abs(vfn-v_guess)));
        v_guess=vfn;
    end
    
    pol_idx=permute(pol_idx,[3 1 2]);
    g=a(pol_idx);
   
    %Setting up initial distribution
    Mu=ones(size(g))/numel(g);
    
    mu_tol=1;
    while mu_tol>1e-08;
        [emp_ind, a_ind ,mass]=find(Mu);
        MuNew=zeros(size(Mu));
        for ii=1:length(emp_ind);
            apr_ind=pol_idx(emp_ind(ii),a_ind(ii));
            MuNew(:,apr_ind)=MuNew(:,apr_ind)+(logzprob(emp_ind(ii),:)*mass(ii))';
        end;
        mu_tol=max(abs(MuNew(:)-Mu(:)));
        
        Mu=MuNew;
    end
    
    aggsav=sum(g(:).*Mu(:));
    
    if aggsav>0;
        kmax=k_guess;
    else
        kmin=k_guess;
    end
display (['k=',num2str(k_guess)])
display (['aggsav=',num2str(aggsav)])
display (['New amax=',num2str(amax)])
display (['new k=',num2str((kmin+kmax)/2)])
display (['vtol=',num2str(v_tol)])
    


  
 
end;
    
  