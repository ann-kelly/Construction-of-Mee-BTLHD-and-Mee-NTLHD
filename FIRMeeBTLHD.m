function [D,ctheta,cnugget]=FIRMeeBTLHD(n,d,theta,nugget,MmLHD)
le1=length(theta);le2=length(nugget);
le=le1*le2;
init=MmLHD;
bestD=zeros(n,d,le);
ropt=zeros(le,d);
for i=1:le
    p1=ceil(i/le2); p2=i-(p1-1)*le2;
    objective =@(r) btFIREEC(init,r,theta(p1),nugget(p2),n,d);
    options = optimoptions('patternsearch','TolFun',10^-6);
    ropt(i,:)=patternsearch(objective,0.5*ones(1,d),[],[],[],[],0.01*ones(1,d),ones(1,d),[],options);
    bestD(:,:,i) = transform(init,ropt(i,:));
end
[D,pos]=caleffi(bestD,theta,nugget,n,d);
c1=ceil(pos/le2); c2=pos-(c1-1)*le2;
ctheta=theta(c1);
cnugget=nugget(c2);


function F_X = btFIREEC(X,r,theta,nugget,n,d)
X=transform(X,r);
storeR=zeros(n,n,d);
for j=1:d 
    storeR(:,:,j) = exp(-(repmat(X(:,j),1,n)-repmat(X(:,j)',n,1)).^2);
end
store=zeros(d,1);
F_X=0;
for j=1:d
    R = prod(storeR(:,:,1:j),3).^(1/(theta*j));
if j<=(d-1)
    store(j)=sum(log(svd(R)+nugget)); 
else
    store(j)=sum(log(svd(R)));
end
    F_X=F_X-store(j);
end
F_X=F_X/d;
end

function TAX=transform(X,r)
[num,dim]=size(X);
TAX=zeros(num,dim);
for j=1:dim
    for k=1:num
        TAX(k,j)=betainv(X(k,j),r(j),r(j));
    end
end
end

function [D,pos]=caleffi(input_designs,theta,nugget,n,d)
ntheta = length(theta); nnugget=length(nugget);
npara = ntheta*nnugget;
fvals = zeros(npara,npara);
for k = 1:npara
    for t = 1:npara
        k1=ceil(k/nnugget); k2=k-(k1-1)*nnugget;
        fvals(t,k) = FIREEC(input_designs(:,:,t),theta(k1),nugget(k2),n,d)-(1-1/d)*n*log(1+nugget(k2));
    end
end
max_at_theta = max(fvals);
fvals = meshgrid(max_at_theta,1:npara)./fvals;
min_at_D = min(fvals,[],2);
[~,pos]=max(min_at_D);
D=input_designs(:,:,pos);  
end

function F_X =  FIREEC(X,theta,nugget,n,d)
storeR=zeros(n,n,d);
for j=1:d 
    storeR(:,:,j) = exp(-(repmat(X(:,j),1,n)-repmat(X(:,j)',n,1)).^2);
end
store=zeros(d,1);
F_X=0;
for j=1:d
    R = prod(storeR(:,:,1:j),3).^(1/(theta*j));
if j<=(d-1)
    store(j)=sum(log(svd(R)+nugget)); 
else
    store(j)=sum(log(svd(R)));
end
    F_X=F_X+store(j);
end
F_X=F_X/d;
end

end