clc;close all;
M=8;
N=4;
HN=1000;
RHN=5*HN;
H = zeros(M,N,RHN);
A=importdata('Hdata.txt'); %%1000行32列
for j=0:(RHN/HN-1)
    for i=1:HN
        H(:,:,i+j*HN)=(reshape(A(i,:),N,M))';%H(:,:,i)代表第i个随机矩阵
    end
end
dnP=2;
nP=-6:dnP:2;
iterN_vector=[1,4,8,12,12];
errnmat=zeros(length(iterN_vector),length(nP));
s=rand(N,RHN);
s(s<0.5)=-1;
s(s>=0.5)=1;
for ii=1:length(iterN_vector)
    for snrn=1:length(nP)
        y=ones(N,RHN);
        r=zeros(M,RHN);
        sigma=sqrt(1/(2*10^(nP(snrn)/10)));
        for i=1:RHN
            r(:,i)=H(:,:,i)*s(:,i);
            r(:,i)=r(:,i)+sigma*randn(M,1);
            if ii==1
                y(:,i)=(H(:,:,i)'*H(:,:,i))^(-1)*H(:,:,i)'*r(:,i);
            end
            if ii>1&&ii<5
                b=H(:,:,i)'*r(:,i);
                A=H(:,:,i)'*H(:,:,i);
                for iterN=1:iterN_vector(ii)
                    d=b-A*y(:,i);
                    w=d'*A*d;
                    t=d'*d/w;
                    y(:,i)=y(:,i)+t*d;
                end
            end
            if ii>=5
                b=H(:,:,i)'*r(:,i);
                A=H(:,:,i)'*H(:,:,i);
                for iterN=1:iterN_vector(ii)
                    d=b-A*y(:,i);
                    w=d'*A*d;
                    x=10^(-5);
                    for iii=1:15
                        x=x*(2-x*w);
                    end
                    t=d'*d*x;
                    y(:,i)=y(:,i)+t*d;
                end
            end
        end
        sr=y;
        sr(sr>0)=1;
        sr(sr<0)=-1;
        errmat=s-sr;
        errn=sum(sum(errmat~=0));
        errnmat(ii,snrn)=errn;
    end
end
BERmat=errnmat./(RHN*N);
figure(1)
%semilogy(nP,BERmat(1,:),'--sk');
%hold on
semilogy(nP,BERmat(2,:),'-*r');
hold on
semilogy(nP,BERmat(3,:),'-ob');
hold on
semilogy(nP,BERmat(4,:),'-^k');
hold on
semilogy(nP,BERmat(5,:),'-pg');
xlabel('Eb/No(dB)');
ylabel('BER');
legend('SD approach with 4 iteration','SD approach with 8 Iteration','SD approach with 12 Iteration','Newton approach with 12 Iteration');


