function D=MeeBTLHD_matern(n,d,theta,nugget,MmLHD)
le1=length(theta);le2=length(nugget);
le=le1*le2;
init=MmLHD;
bestD=zeros(n,d,le);
ropt=zeros(1,le);
E=logical(fullfact(2*ones(1,d))-1);
E(1,:)=[];
N=size(E,1);
theta=theta.*0.135;
for i=1:le
    p1=ceil(i/le2); p2=i-(p1-1)*le2;
    objective = @(r) btEEC(init,r,theta(p1),nugget(p2),n,d,E,N);
    options = optimoptions('patternsearch','TolFun',10^-6);
    ropt(i)=patternsearch(objective,0.5,[],[],[],[],0.01,1,[],options);
    bestD(:,:,i) = transform(init,ropt(i));
end
D=caleffi(bestD,theta,nugget,n,d,E,N);

function F_X = btEEC(X,r,theta,nugget,n,d,E,N)
X=transform(X,r);
storeR=zeros(n,n,d);
for j=1:d 
    storeR(:,:,j) = (repmat(X(:,j),1,n)-repmat(X(:,j)',n,1)).^2;
end
store=zeros(N,1);
F_X=0;
for j=1:N
    s=sum(E(j,:));
    Distance=sqrt(sum(storeR(:,:,E(j,:)),3))./sqrt(theta*s);
    R = (1+Distance+Distance.^2/3).*exp(-Distance);
if s<=(d-1)
    store(j)=sum(log(svd(R)+nugget)); 
else
    store(j)=sum(log(svd(R)));
end
    F_X=F_X-store(j);
end
F_X=F_X/N;
end

function TAX = transform(X,r)
    TAX=betainv(X,r,r);
end


function D=caleffi(input_designs,theta,nugget,n,d,E,N)
ntheta = length(theta); nnugget=length(nugget);
npara = ntheta*nnugget;
fvals = zeros(npara,npara);
for k = 1:npara
    for t = 1:npara
        k1=ceil(k/nnugget); k2=k-(k1-1)*nnugget;
        fvals(t,k) = EEC(input_designs(:,:,t),theta(k1),nugget(k2),n,d,E,N)-(1-1/N)*n*log(1+nugget(k2));
    end
end
max_at_theta = max(fvals);
fvals = meshgrid(max_at_theta,1:npara)./fvals;
min_at_D = min(fvals,[],2);
[~,pos]=max(min_at_D);
D=input_designs(:,:,pos);  
end



function F_X = EEC(x,theta,nugget,n,d,E,N)
X = reshape(x,[n,d]);
storeR=zeros(n,n,d);
for j=1:d 
    storeR(:,:,j) = (repmat(X(:,j),1,n)-repmat(X(:,j)',n,1)).^2;
end
store=zeros(N,1);
F_X=0;
for j=1:N
    s=sum(E(j,:));
    Distance=sqrt(sum(storeR(:,:,E(j,:)),3))./sqrt(theta*s);
    R = (1+Distance+Distance.^2/3).*exp(-Distance);
if s<=(d-1)
    store(j)=sum(log(svd(R)+nugget)); 
else
    store(j)=sum(log(svd(R))); 
end
    F_X=F_X+store(j);
end
F_X=F_X/N;
end

end