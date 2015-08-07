function mask2=post_processing(mask1_1)
% L_1= bwlabel(mask1_1);
% max_region_1=max(max(L_1));
% if max_region_1>10000
%     si=[];
%     for i=1:max_region_1
%         si(i)=size(find(L_1==i),1);
%     end
%     [max_s,index]=sort(si,'descend');
%     T=max_s(1)*0.5;
%     for i=1:max_region_1
%         if si(i)<=T
%             mask1_1(find(L_1==i))=0;
%         end
%     end
% end

L_1= bwlabel(mask1_1);
max_region_1=max(max(L_1));
if max_region_1~=1 
    si=[];
    for i=1:max_region_1
        si(i)=size(find(L_1==i),1);
    end
    [max_s,index]=sort(si,'descend');
    T=max_s(1)*0.2;
    for i=1:max_region_1
        if si(i)<=T
            mask1_1(find(L_1==i))=0;
        end
    end
end
mask2=mask1_1;