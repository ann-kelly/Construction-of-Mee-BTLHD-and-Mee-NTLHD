[n,d]=size(MmLHD);
Combinations=fullfact(2*ones(1,d))-1;
Combinations=logical(Combinations(sum(Combinations,2)>0,:));
m=size(Combinations,1);
plot_theta=0.05:0.01:2;
plot_n=length(plot_theta);
nugget=10^-3;
LogDet1=zeros(m,plot_n);LogDet2=zeros(m,plot_n);
LogDet3=zeros(m,plot_n);LogDet4=zeros(m,plot_n);
LogDet5=zeros(m,plot_n);LogDet6=zeros(m,plot_n);
for i=1:m
for j=1:plot_n
    k=sum(Combinations(i,:));
    if k==d
        R1=LG(MeeNTLHDm,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval1=svd(R1);LogDet1(i,j)=sum(log(Eigval1));
        R2=LG(MeeBTLHDm,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval2=svd(R2);LogDet2(i,j)=sum(log(Eigval2));
        R3=LG(uniform,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval3=svd(R3);LogDet3(i,j)=sum(log(Eigval3));
        R4=LG(MmLHD,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval4=svd(R4);LogDet4(i,j)=sum(log(Eigval4));
        R5=LG(MaxPro,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval5=svd(R5);LogDet5(i,j)=sum(log(Eigval5));
        R6=LG(GMLHD,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval6=svd(R6);LogDet6(i,j)=sum(log(Eigval6));
    else
        R1=LG(MeeNTLHDm,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval1=svd(R1)+nugget;LogDet1(i,j)=sum(log(Eigval1));
        R2=LG(MeeBTLHDm,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval2=svd(R2)+nugget;LogDet2(i,j)=sum(log(Eigval2));
        R3=LG(uniform,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval3=svd(R3)+nugget;LogDet3(i,j)=sum(log(Eigval3));
        R4=LG(MmLHD,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval4=svd(R4)+nugget;LogDet4(i,j)=sum(log(Eigval4));
        R5=LG(MaxPro,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval5=svd(R5)+nugget;LogDet5(i,j)=sum(log(Eigval5));
        R6=LG(GMLHD,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval6=svd(R6)+nugget;LogDet6(i,j)=sum(log(Eigval6));
    end   
end
end
a=zeros(1,d);
averent1=zeros(d,plot_n);averent2=zeros(d,plot_n);
averent3=zeros(d,plot_n);averent4=zeros(d,plot_n);
averent5=zeros(d,plot_n);averent6=zeros(d,plot_n);
for i=1:m
    j=sum(Combinations(i,:));
    a(j)=a(j)+1;
    averent1(j,:)= averent1(j,:)+LogDet1(i,:);
    averent2(j,:)= averent2(j,:)+LogDet2(i,:);
    averent3(j,:)= averent3(j,:)+LogDet3(i,:);
    averent4(j,:)= averent4(j,:)+LogDet4(i,:);
    averent5(j,:)= averent5(j,:)+LogDet5(i,:);
    averent6(j,:)= averent6(j,:)+LogDet6(i,:);
end
for i=1:(d-1)
    figure(i);
    ax1=subplot(2,1,1);
    set(gca,'XTick',[0 0.1 0.2:0.2:2]);
    hold on
    plot(ax1,plot_theta,averent1(i,:)/a(i)-averent3(i,:)/a(i),'-sk','MarkerIndices',1:10:200);
    plot(ax1,plot_theta,averent1(i,:)/a(i)-averent4(i,:)/a(i),'--+k','MarkerIndices',1:10:200);
    plot(ax1,plot_theta,averent1(i,:)/a(i)-averent5(i,:)/a(i),':*k','MarkerIndices',1:10:200);
    plot(ax1,plot_theta,averent1(i,:)/a(i)-averent6(i,:)/a(i),'-.xk','MarkerIndices',1:10:200);
    ylabel('Relative average entropy');
    title([num2str(i),'D projection'],'FontWeight','bold');
end
figure(1);legend('uniform','MmLHD','MaxPro','GMLHD');
figure(d+1);
ax1=subplot(2,1,1);
set(gca,'XTick',[0 0.1 0.2:0.2:2]);
hold on
plot(ax1,plot_theta,averent1(1,:)/a(1)-averent2(1,:)/a(1),'-ok','MarkerIndices',1:10:200);
plot(ax1,plot_theta,averent1(2,:)/a(2)-averent2(2,:)/a(2),'--+k','MarkerIndices',1:10:200);
plot(ax1,plot_theta,averent1(3,:)/a(3)-averent2(3,:)/a(3),'-.xk','MarkerIndices',1:10:200);
plot(ax1,plot_theta,averent1(4,:)/a(4)-averent2(4,:)/a(4),':*k','MarkerIndices',1:10:200);
plot(ax1,plot_theta,averent1(5,:)/a(5)-averent2(5,:)/a(5),'-sk','MarkerIndices',1:10:200);
ylabel('Relative average entropy');
title('Relative average entropy for Mee-BTLHD','FontWeight','bold');


nugget=10^-1;
for i=1:m
for j=1:plot_n
    k=sum(Combinations(i,:));
    if k<d
        R1=LG(MeeNTLHDm,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval1=svd(R1)+nugget;LogDet1(i,j)=sum(log(Eigval1));
        R2=LG(MeeBTLHDm,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval2=svd(R2)+nugget;LogDet2(i,j)=sum(log(Eigval2));
        R3=LG(uniform,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval3=svd(R3)+nugget;LogDet3(i,j)=sum(log(Eigval3));
        R4=LG(MmLHD,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval4=svd(R4)+nugget;LogDet4(i,j)=sum(log(Eigval4));
        R5=LG(MaxPro,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval5=svd(R5)+nugget;LogDet5(i,j)=sum(log(Eigval5));
        R6=LG(GMLHD,n,k,Combinations(i,:),plot_theta(j)*0.135);Eigval6=svd(R6)+nugget;LogDet6(i,j)=sum(log(Eigval6));
    end   
end
end
a=zeros(1,d);
averent1=zeros(d,plot_n);averent2=zeros(d,plot_n);
averent3=zeros(d,plot_n);averent4=zeros(d,plot_n);
averent5=zeros(d,plot_n);averent6=zeros(d,plot_n);
for i=1:m
    j=sum(Combinations(i,:));
    a(j)=a(j)+1;
    averent1(j,:)= averent1(j,:)+LogDet1(i,:);
    averent2(j,:)= averent2(j,:)+LogDet2(i,:);
    averent3(j,:)= averent3(j,:)+LogDet3(i,:);
    averent4(j,:)= averent4(j,:)+LogDet4(i,:);
    averent5(j,:)= averent5(j,:)+LogDet5(i,:);
    averent6(j,:)= averent6(j,:)+LogDet6(i,:);
end
for i=1:(d-1)
    figure(i);
    ax2=subplot(2,1,2);
    set(gca,'XTick',[0 0.1 0.2:0.2:2]);
    hold on
    plot(ax2,plot_theta,averent1(i,:)/a(i)-averent3(i,:)/a(i),'-sk','MarkerIndices',1:10:200);
    plot(ax2,plot_theta,averent1(i,:)/a(i)-averent4(i,:)/a(i),'--+k','MarkerIndices',1:10:200);
    plot(ax2,plot_theta,averent1(i,:)/a(i)-averent5(i,:)/a(i),':*k','MarkerIndices',1:10:200);
    plot(ax2,plot_theta,averent1(i,:)/a(i)-averent6(i,:)/a(i),'-.xk','MarkerIndices',1:10:200);
    xlabel('\theta');
    ylabel('Relative average entropy');
end
figure(d);
set(gca,'XTick',[0 0.1 0.2:0.2:2]);
hold on
plot(plot_theta,averent1(d,:)/a(d)-averent3(d,:)/a(d),'-sk','MarkerIndices',1:10:200);
plot(plot_theta,averent1(d,:)/a(d)-averent4(d,:)/a(d),'--+k','MarkerIndices',1:10:200);
plot(plot_theta,averent1(d,:)/a(d)-averent5(d,:)/a(d),':*k','MarkerIndices',1:10:200);
plot(plot_theta,averent1(d,:)/a(d)-averent6(d,:)/a(d),'-.xk','MarkerIndices',1:10:200);
xlabel('\theta');
ylabel('Relative average entropy');
title('Unprojected','FontWeight','bold');
figure(d+1);
ax2=subplot(2,1,2);
set(gca,'XTick',[0 0.1 0.2:0.2:2]);
hold on
plot(ax2,plot_theta,averent1(1,:)/a(1)-averent2(1,:)/a(1),'-ok','MarkerIndices',1:10:200);
plot(ax2,plot_theta,averent1(2,:)/a(2)-averent2(2,:)/a(2),'--+k','MarkerIndices',1:10:200);
plot(ax2,plot_theta,averent1(3,:)/a(3)-averent2(3,:)/a(3),'-.xk','MarkerIndices',1:10:200);
plot(ax2,plot_theta,averent1(4,:)/a(4)-averent2(4,:)/a(4),':*k','MarkerIndices',1:10:200);
plot(ax2,plot_theta,averent1(5,:)/a(5)-averent2(5,:)/a(5),'-sk','MarkerIndices',1:10:200);
legend('1D','2D','3D','4D','5D');
xlabel('\theta');
ylabel('Relative average entropy');


function R=LG(D,n,k,comb,theta)
pD=D(:,comb);Distance=pdist2(pD,pD)./sqrt(theta*k);
R=(1+Distance+Distance.^2/3).*exp(-Distance);
end