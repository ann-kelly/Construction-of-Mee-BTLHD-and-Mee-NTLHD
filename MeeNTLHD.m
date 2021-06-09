function D=MeeNTLHD(n,d,theta,nugget,MmLHD)
le1=length(theta);le2=length(nugget);
le=le1*le2;
m=floor(n/2);
init=MmLHD;
bestD=zeros(n,d,le);
ropt=zeros(le,m);
E=logical(fullfact(2*ones(1,d))-1);
E(1,:)=[];
N=size(E,1);
for i=1:le
    p1=ceil(i/le2); p2=i-(p1-1)*le2;
    objective =@(r) ntEEC(init,r,theta(p1),nugget(p2),n,d,E,N);
    A=eye(m);
    A(repmat((1:m)',1,m)==repmat((1:m)-1,m,1))=-1;
    b=1/n*ones(m,1);
    b(m)=0.5/n;
    b=b-10^-3;
    options = optimoptions('patternsearch','TolFun',10^-6,'MaxFunEvals',2000*m);
    ropt(i,:)=patternsearch(objective,zeros(1,m),A,b,[],[],zeros(1,m),[],[],options);
    bestD(:,:,i) = transform(init,ropt(i,:));
end
D=caleffi(bestD,theta,nugget,n,d,E,N);

function F_X = ntEEC(X,r,theta,nugget,n,d,E,N)
X=transform(X,r);
storeR=zeros(n,n,d);
for j=1:d 
    storeR(:,:,j) = exp(-(repmat(X(:,j),1,n)-repmat(X(:,j)',n,1)).^2);
end
store=zeros(N,1);
F_X=0;
for j=1:N
    s=sum(E(j,:));
    R = prod(storeR(:,:,E(j,:)),3).^(1/(theta*s));
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
    num=size(X,1);
    TAX=X;
    TAX2=ceil(X*num);
    for j=1:floor(num/2)
        I=(TAX2-j)==0;
        TAX(I)=TAX(I)-r(end+1-j);
        I2=(TAX2-(num+1-j))==0;
        TAX(I2)=TAX(I2)+r(end+1-j);
    end
end

function D=caleffi(input_designs,theta,nugget,n,d,E,N)
ntheta = length(theta);nnugget=length(nugget);
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
    storeR(:,:,j) = exp(-(repmat(X(:,j),1,n)-repmat(X(:,j)',n,1)).^2);
end
store=zeros(N,1);
F_X=0;
for j=1:N
    s=sum(E(j,:));
    R = prod(storeR(:,:,E(j,:)),3).^(1/(theta*s));
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