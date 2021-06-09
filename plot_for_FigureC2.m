[n,d]=size(MmLHD);
order1=orderfactor(GMLHD,ctheta,cnugget,d);
GMLHD_new=GMLHD(:,order1);
order2=orderfactor(MaxPro,ctheta,cnugget,d);
MaxPro_new=MaxPro(:,order2);
order3=orderfactor(MmLHD,ctheta,cnugget,d);
MmLHD_new=MmLHD(:,order3);
order4=orderfactor(uniform,ctheta,cnugget,d);
uniform_new=uniform(:,order4);
plot_theta=0.05:0.01:2;
plot_n=length(plot_theta);
nugget=10^-3;
LogDet1=zeros(d,plot_n);
LogDet2=zeros(d,plot_n);LogDet3=zeros(d,plot_n);
LogDet4=zeros(d,plot_n);LogDet5=zeros(d,plot_n);
for i=1:d
for j=1:plot_n
    if i==d
        D1=FIRMeeBTLHDg(:,1:i);Distances1=pdist2(D1,D1);R1=exp(-Distances1.^2/(plot_theta(j)*i));Eigval1=svd(R1);LogDet1(i,j)=sum(log(Eigval1));
        D2=uniform_new(:,1:i);Distances2=pdist2(D2,D2);R2=exp(-Distances2.^2/(plot_theta(j)*i));Eigval2=svd(R2);LogDet2(i,j)=sum(log(Eigval2));
        D3=MmLHD_new(:,1:i);Distances3=pdist2(D3,D3);R3=exp(-Distances3.^2/(plot_theta(j)*i));Eigval3=svd(R3);LogDet3(i,j)=sum(log(Eigval3));
        D4=MaxPro_new(:,1:i);Distances4=pdist2(D4,D4);R4=exp(-Distances4.^2/(plot_theta(j)*i));Eigval4=svd(R4);LogDet4(i,j)=sum(log(Eigval4));
        D5=GMLHD_new(:,1:i);Distances5=pdist2(D5,D5);R5=exp(-Distances5.^2/(plot_theta(j)*i));Eigval5=svd(R5);LogDet5(i,j)=sum(log(Eigval5));
    else
        D1=FIRMeeBTLHDg(:,1:i);Distances1=pdist2(D1,D1);R1=exp(-Distances1.^2/(plot_theta(j)*i));Eigval1=svd(R1)+nugget;LogDet1(i,j)=sum(log(Eigval1));
        D2=uniform_new(:,1:i);Distances2=pdist2(D2,D2);R2=exp(-Distances2.^2/(plot_theta(j)*i));Eigval2=svd(R2)+nugget;LogDet2(i,j)=sum(log(Eigval2));
        D3=MmLHD_new(:,1:i);Distances3=pdist2(D3,D3);R3=exp(-Distances3.^2/(plot_theta(j)*i));Eigval3=svd(R3)+nugget;LogDet3(i,j)=sum(log(Eigval3));
        D4=MaxPro_new(:,1:i);Distances4=pdist2(D4,D4);R4=exp(-Distances4.^2/(plot_theta(j)*i));Eigval4=svd(R4)+nugget;LogDet4(i,j)=sum(log(Eigval4));
        D5=GMLHD_new(:,1:i);Distances5=pdist2(D5,D5);R5=exp(-Distances5.^2/(plot_theta(j)*i));Eigval5=svd(R5)+nugget;LogDet5(i,j)=sum(log(Eigval5));
    end   
end
end
for i=1:(d-1)
    figure(i);
    ax1=subplot(2,1,1);
    set(gca,'XTick',[0 0.1 0.2:0.2:2]);
    hold on
    plot(ax1,plot_theta,LogDet1(i,:)-LogDet2(i,:),'-sk','MarkerIndices',1:10:200);
    plot(ax1,plot_theta,LogDet1(i,:)-LogDet3(i,:),'--+k','MarkerIndices',1:10:200);
    plot(ax1,plot_theta,LogDet1(i,:)-LogDet4(i,:),':*k','MarkerIndices',1:10:200);
    plot(ax1,plot_theta,LogDet1(i,:)-LogDet5(i,:),'-.xk','MarkerIndices',1:10:200);
    ylabel('Relative average entropy');
    title([num2str(i),'D projection'],'FontWeight','bold');
end
figure(1);legend('uniform','MmLHD','MaxPro','GMLHD');

nugget=10^-1;
LogDet1=zeros(d,plot_n);
LogDet2=zeros(d,plot_n);LogDet3=zeros(d,plot_n);
LogDet4=zeros(d,plot_n);LogDet5=zeros(d,plot_n);
for i=1:d
for j=1:plot_n
    if i==d
        D1=FIRMeeBTLHDg(:,1:i);Distances1=pdist2(D1,D1);R1=exp(-Distances1.^2/(plot_theta(j)*i));Eigval1=svd(R1);LogDet1(i,j)=sum(log(Eigval1));
        D2=uniform_new(:,1:i);Distances2=pdist2(D2,D2);R2=exp(-Distances2.^2/(plot_theta(j)*i));Eigval2=svd(R2);LogDet2(i,j)=sum(log(Eigval2));
        D3=MmLHD_new(:,1:i);Distances3=pdist2(D3,D3);R3=exp(-Distances3.^2/(plot_theta(j)*i));Eigval3=svd(R3);LogDet3(i,j)=sum(log(Eigval3));
        D4=MaxPro_new(:,1:i);Distances4=pdist2(D4,D4);R4=exp(-Distances4.^2/(plot_theta(j)*i));Eigval4=svd(R4);LogDet4(i,j)=sum(log(Eigval4));
        D5=GMLHD_new(:,1:i);Distances5=pdist2(D5,D5);R5=exp(-Distances5.^2/(plot_theta(j)*i));Eigval5=svd(R5);LogDet5(i,j)=sum(log(Eigval5));
    else
        D1=FIRMeeBTLHDg(:,1:i);Distances1=pdist2(D1,D1);R1=exp(-Distances1.^2/(plot_theta(j)*i));Eigval1=svd(R1)+nugget;LogDet1(i,j)=sum(log(Eigval1));
        D2=uniform_new(:,1:i);Distances2=pdist2(D2,D2);R2=exp(-Distances2.^2/(plot_theta(j)*i));Eigval2=svd(R2)+nugget;LogDet2(i,j)=sum(log(Eigval2));
        D3=MmLHD_new(:,1:i);Distances3=pdist2(D3,D3);R3=exp(-Distances3.^2/(plot_theta(j)*i));Eigval3=svd(R3)+nugget;LogDet3(i,j)=sum(log(Eigval3));
        D4=MaxPro_new(:,1:i);Distances4=pdist2(D4,D4);R4=exp(-Distances4.^2/(plot_theta(j)*i));Eigval4=svd(R4)+nugget;LogDet4(i,j)=sum(log(Eigval4));
        D5=GMLHD_new(:,1:i);Distances5=pdist2(D5,D5);R5=exp(-Distances5.^2/(plot_theta(j)*i));Eigval5=svd(R5)+nugget;LogDet5(i,j)=sum(log(Eigval5));
    end   
end
end
for i=1:(d-1)
    figure(i);
    ax2=subplot(2,1,2);
    set(gca,'XTick',[0 0.1 0.2:0.2:2]);
    hold on
    plot(ax2,plot_theta,LogDet1(i,:)-LogDet2(i,:),'-sk','MarkerIndices',1:10:200);
    plot(ax2,plot_theta,LogDet1(i,:)-LogDet3(i,:),'--+k','MarkerIndices',1:10:200);
    plot(ax2,plot_theta,LogDet1(i,:)-LogDet4(i,:),':*k','MarkerIndices',1:10:200);
    plot(ax2,plot_theta,LogDet1(i,:)-LogDet5(i,:),'-.xk','MarkerIndices',1:10:200);
    xlabel('\theta');
    ylabel('Relative average entropy');
end
figure(d);
set(gca,'XTick',[0 0.1 0.2:0.2:2]);
hold on
plot(plot_theta,LogDet1(d,:)-LogDet2(d,:),'-sk','MarkerIndices',1:10:200);
plot(plot_theta,LogDet1(d,:)-LogDet3(d,:),'--+k','MarkerIndices',1:10:200);
plot(plot_theta,LogDet1(d,:)-LogDet4(d,:),':*k','MarkerIndices',1:10:200);
plot(plot_theta,LogDet1(d,:)-LogDet5(d,:),'-.xk','MarkerIndices',1:10:200);
xlabel('\theta');
ylabel('Relative average entropy');
title('Unprojected','FontWeight','bold');

function choose_d_result = orderfactor(x,theta,nugget,d)
choose_d = zeros(1, d);
det_a = zeros(d, 1);
for i = 1 : d
    distances=pdist2(x(:,i),x(:,i));
    matrix_a=exp(-distances.^2/theta);
    det_a(i)=sum(log(svd(matrix_a)+nugget));
end
[~,d_max] = max(det_a);
choose_d(1) = d_max;
for i = 2 : d
   if i<=(d-1)
    det_new = -inf(d, 1);
    selected_d = choose_d(1 : i);
    for k = 1 : d
        if any(k == choose_d(1 : i - 1))
            continue;
        end
        selected_d(i) = k;
        distances=pdist2(x(:,selected_d),x(:,selected_d));
        matrix_b=exp(-distances.^2/(theta*i));
        det_new(k)=sum(log(svd(matrix_b)+nugget));
        
    end
    [~,d_max] = max(det_new);
    choose_d(i) = d_max;
   else
    det_new = -inf(d, 1);
    selected_d = choose_d(1 : i);
    for k = 1 : d
        if any(k == choose_d(1 : i - 1))
            continue;
        end
        selected_d(i) = k;
        distances=pdist2(x(:,selected_d),x(:,selected_d));
        matrix_b=exp(-distances.^2/(theta*i));
        det_new(k)=sum(log(svd(matrix_b)));
        
    end
    [~,d_max] = max(det_new);
    choose_d(i) = d_max;
   end
end
choose_d_result = choose_d;
end