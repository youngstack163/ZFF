clc
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
dnP=0.1;
nP=-6:dnP:2;
iterN_vector=[12,1,4,8,12];
errnmat=zeros(length(iterN_vector),length(nP));
s=rand(N,RHN);
s(s<0.5)=-1;
s(s>=0.5)=1;
for ii=1:length(iterN_vector)
    for snrn=1:length(nP)
        y=ones(N,RHN);
        r=zeros(M,RHN);
        for i=1:RHN
            r(:,i)=H(:,:,i)*s(:,i);
            r(:,i)=awgn(r(:,i),nP(snrn));
            if ii==1
                y(:,i)=s(:,i);
            end
            if ii==2
                y(:,i)=(H(:,:,i)'*H(:,:,i))^(-1)*H(:,:,i)'*r(:,i);
            end
            if ii>2
                b=H(:,:,i)'*r(:,i);
                A=H(:,:,i)'*H(:,:,i);
                for iterN=1:iterN_vector(ii)
                    d=b-A*y(:,i);
                    w=d'*A*d;
                    t=d'*d/w;
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
plot(nP,BERmat(1,:),'-+');
hold on
plot(nP,BERmat(2,:),'-*');
hold on
plot(nP,BERmat(3,:),'-x');
hold on
plot(nP,BERmat(4,:),'-d');
hold on
plot(nP,BERmat(5,:),'-^');
grid on

