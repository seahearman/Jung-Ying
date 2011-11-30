clc;clear all;
H=zeros(8,5);
H(1,:)=[0 0 0 0 0];
H(2,:)=[0 0 0 0 1];
H(3,:)=[0 1 0 0 0];
H(4,:)=[0 1 0 0 1];
H(5,:)=[0 1 1 0 0];
H(6,:)=[0 1 1 1 1];
H(7,:)=[1 0 0 0 0];
H(8,:)=[1 0 0 0 1];
R=eye(8);

for i=1:8
    for j=i+1:8
        R(i,j)=sum(H(i,:)==H(j,:))/5;
        R(j,i)=R(i,j);
    end
end
disp(R);
dlmwrite('R.txt',R)
        
