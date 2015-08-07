function [phi1,flag,ratio]=rescale_phi(phi,phi_ori,super,nn)
phi1=phi;
% for i=1:length(phi)
i=nn;
    m1=phi{i}<0;
    m2=phi_ori{i}<0;
    m=m1-m2;
    ratio=length(find(m~=0))/(size(m1,1)*size(m1,2));
    if ratio<0.001
        flag=0;
%         phi1{i}=phi{i};
    else
        flag=1;
%         idx=find(super(i).vec==-1);
%         b=[];
%         for j=1:length(idx)
%             idx2=find(super(i).label==idx(j));
%             b(j)=sum(m(idx2))>0;
%             m1(idx2)=b(j);
%         end
%         phi1{i}=mask2phi(m1);
    end
% end
            