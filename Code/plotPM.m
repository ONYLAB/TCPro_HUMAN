for i = 1:33
    pM(i,:)=reshape(Y1(i,(12+7*N):(11+13*N)),N,6);
    pM_NUMBER(i)=EpitopeRatio*sum(pM(i,:)*ones(6,1)*1E-12*NA);
end
plot(T1,pM_NUMBER)


for i = 1:33
    pM(i,:)=reshape(Y1(i,(12+7*N):(11+13*N)),N,6);
    pM_NUMBER(i)=sum(pM(i,:)*ones(6,1)*1E-12*NA);
end
plot(T1,pM_NUMBER)