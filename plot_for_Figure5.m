[n,d]=size(MmLHD);
gn=10;
for i=1:d
    eval(['g',num2str(i),'=','(fullfact(gn*ones(1,i))-1)./(gn-1)']);
end
Combinations=fullfact(2*ones(1,d))-1;
Combinations=logical(Combinations(sum(Combinations,2)>0,:));
m=size(Combinations,1);
plot_theta=0.05:0.01:0.8;
plot_n=length(plot_theta);


nugget=10^-3;
LogDet1=zeros(m,plot_n);LogDet2=zeros(m,plot_n);
LogDet3=zeros(m,plot_n);LogDet4=zeros(m,plot_n);
LogDet5=zeros(m,plot_n);LogDet6=zeros(m,plot_n);
for i=1:m
for j=1:plot_n
    k=sum(Combinations(i,:));
    if k==d
        LogDet1(i,j)=povar(MeeNTLHDg,n,k,Combinations(i,:),plot_theta(j),0,eval(['g',num2str(k)]));
        LogDet2(i,j)=povar(MeeBTLHDg,n,k,Combinations(i,:),plot_theta(j),0,eval(['g',num2str(k)]));
        LogDet3(i,j)=povar(uniform,n,k,Combinations(i,:),plot_theta(j),0,eval(['g',num2str(k)]));
        LogDet4(i,j)=povar(MmLHD,n,k,Combinations(i,:),plot_theta(j),0,eval(['g',num2str(k)]));
        LogDet5(i,j)=povar(MaxPro,n,k,Combinations(i,:),plot_theta(j),0,eval(['g',num2str(k)]));
        LogDet6(i,j)=povar(GMLHD,n,k,Combinations(i,:),plot_theta(j),0,eval(['g',num2str(k)]));
    else
        LogDet1(i,j)=povar(MeeNTLHDg,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet2(i,j)=povar(MeeBTLHDg,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet3(i,j)=povar(uniform,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet4(i,j)=povar(MmLHD,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet5(i,j)=povar(MaxPro,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet6(i,j)=povar(GMLHD,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
    end   
end
end
a=zeros(1,d);
averent1=zeros(d,plot_n);averent2=zeros(d,plot_n);
averent3=zeros(d,plot_n);averent4=zeros(d,plot_n);
averent5=zeros(d,plot_n); averent6=zeros(d,plot_n);
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
    axis([0.05 0.8 -inf inf]);
    set(gca,'XTick',[0.05 0.1:0.1:0.8]);
    hold on
    plot(ax1,plot_theta,averent3(i,:)/a(i)-averent1(i,:)/a(i),'-sk','MarkerIndices',1:7:70);
    plot(ax1,plot_theta,averent4(i,:)/a(i)-averent1(i,:)/a(i),'--+k','MarkerIndices',1:7:70);
    plot(ax1,plot_theta,averent5(i,:)/a(i)-averent1(i,:)/a(i),':*k','MarkerIndices',1:7:70);
    plot(ax1,plot_theta,averent6(i,:)/a(i)-averent1(i,:)/a(i),'-.xk','MarkerIndices',1:7:70);
    ylabel('Relative MPV');
    title([num2str(i),'D projection'],'FontWeight','bold');
end
figure(1);legend('uniform','MmLHD','MaxPro','GMLHD');
figure(d+1);
ax1=subplot(2,1,1);
axis([0.05 0.8 -inf inf]);
set(gca,'XTick',[0.05 0.1:0.1:0.8]);
hold on
plot(ax1,plot_theta,averent2(1,:)/a(1)-averent1(1,:)/a(1),'-ok','MarkerIndices',1:7:70);
plot(ax1,plot_theta,averent2(2,:)/a(2)-averent1(2,:)/a(2),'--+k','MarkerIndices',1:7:70);
plot(ax1,plot_theta,averent2(3,:)/a(3)-averent1(3,:)/a(3),'-.xk','MarkerIndices',1:7:70);
plot(ax1,plot_theta,averent2(4,:)/a(4)-averent1(4,:)/a(4),':*k','MarkerIndices',1:7:70);
plot(ax1,plot_theta,averent2(5,:)/a(5)-averent1(5,:)/a(5),'-sk','MarkerIndices',1:7:70);
ylabel('Relative MPV');
title('Relative MPV for Mee-BTLHD','FontWeight','bold');



nugget=10^-1;
for i=1:m
for j=1:plot_n
    k=sum(Combinations(i,:));
    if k<d
        LogDet1(i,j)=povar(MeeNTLHDg,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet2(i,j)=povar(MeeBTLHDg,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet3(i,j)=povar(uniform,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet4(i,j)=povar(MmLHD,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet5(i,j)=povar(MaxPro,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
        LogDet6(i,j)=povar(GMLHD,n,k,Combinations(i,:),plot_theta(j),nugget,eval(['g',num2str(k)]));
    end   
end
end
a=zeros(1,d);
averent1=zeros(d,plot_n);averent2=zeros(d,plot_n);
averent3=zeros(d,plot_n);averent4=zeros(d,plot_n);
averent5=zeros(d,plot_n); averent6=zeros(d,plot_n);
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
    axis([0.05 0.8 -inf inf]);
    set(gca,'XTick',[0.05 0.1:0.1:0.8]);
    hold on
    plot(ax2,plot_theta,averent3(i,:)/a(i)-averent1(i,:)/a(i),'-sk','MarkerIndices',1:7:70);
    plot(ax2,plot_theta,averent4(i,:)/a(i)-averent1(i,:)/a(i),'--+k','MarkerIndices',1:7:70);
    plot(ax2,plot_theta,averent5(i,:)/a(i)-averent1(i,:)/a(i),':*k','MarkerIndices',1:7:70);
    plot(ax2,plot_theta,averent6(i,:)/a(i)-averent1(i,:)/a(i),'-.xk','MarkerIndices',1:7:70);
    xlabel('\theta');
    ylabel('Relative MPV');
end
figure(d);
axis([0.05 0.8 -inf inf]);
set(gca,'XTick',[0.05 0.1:0.1:0.8]);
hold on
plot(plot_theta,averent3(d,:)/a(d)-averent1(d,:)/a(d),'-sk','MarkerIndices',1:7:70);
plot(plot_theta,averent4(d,:)/a(d)-averent1(d,:)/a(d),'--+k','MarkerIndices',1:7:70);
plot(plot_theta,averent5(d,:)/a(d)-averent1(d,:)/a(d),':*k','MarkerIndices',1:7:70);
plot(plot_theta,averent6(d,:)/a(d)-averent1(d,:)/a(d),'-.xk','MarkerIndices',1:7:70);
xlabel('\theta');
ylabel('Relative MPV');
title('Unprojected','FontWeight','bold');
figure(d+1);
ax2=subplot(2,1,2);
axis([0.05 0.8 -inf inf]);
set(gca,'XTick',[0.05 0.1:0.1:0.8]);
hold on
plot(ax2,plot_theta,averent2(1,:)/a(1)-averent1(1,:)/a(1),'-ok','MarkerIndices',1:7:70);
plot(ax2,plot_theta,averent2(2,:)/a(2)-averent1(2,:)/a(2),'--+k','MarkerIndices',1:7:70);
plot(ax2,plot_theta,averent2(3,:)/a(3)-averent1(3,:)/a(3),'-.xk','MarkerIndices',1:7:70);
plot(ax2,plot_theta,averent2(4,:)/a(4)-averent1(4,:)/a(4),':*k','MarkerIndices',1:7:70);
plot(ax2,plot_theta,averent2(5,:)/a(5)-averent1(5,:)/a(5),'-sk','MarkerIndices',1:7:70);
legend('1D','2D','3D','4D','5D');
xlabel('\theta');
ylabel('Relative MPV');




function mpv=povar(D,n,k,comb,theta,nugget,g)
pD=D(:,comb);Distances=pdist2(pD,pD);R=exp(-Distances.^2/(theta*k))+nugget*eye(n);
num=size(g,1);sigma=zeros(1,num);
for i=1:num
    distances=pdist2(g(i,:),pD);
    r=exp(-distances.^2/(theta*k));
    sigma(i)=r/R*r';
end
mpv=max(1+nugget-sigma);
end