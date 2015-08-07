function get_save_results(img_path,phi1,img1,para_energy,img_name)
mkdir([img_path,'results']);
mask1_1=zeros(size(phi1,1),size(phi1,2));
mask1_1(find(phi1<=0))=1;
mask1_1=logical(mask1_1);
imwrite(mask1_1,[img_path,'results_',num2str(para_energy(1)),'_',num2str(para_energy(2)),'/',img_name,'_mask1_ori.bmp'],'bmp');

img1_temp=[];
img1_temp=double(img1);
for j=1:size(img1,1)
    for k=1:size(img1,2)
        if mask1_1(j,k)==0
            img1_temp(j,k,1)=0;
            img1_temp(j,k,2)=0;
            img1_temp(j,k,3)=255;
        end
    end
end
imwrite(uint8(img1_temp),[img_path,'results_',num2str(para_energy(1)),'_',num2str(para_energy(2)),'/',img_name,'_seg1_ori.bmp'],'bmp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%post-processing
L_1= bwlabel(mask1_1);
max_region_1=max(max(L_1));
if max_region_1>100000
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
% [img_paht,'/results_',num2str(para_energy(1)),num2str(para_energy(2))]
imwrite(mask1_1,[img_path,'results_',num2str(para_energy(1)),'_',num2str(para_energy(2)),'/',img_name,'_mask1_post.bmp'],'bmp');
phi11=mask1_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%save the results
img1_temp=[];
img1_temp=img1;
for j=1:size(img1,1)
    for k=1:size(img1,2)
        if mask1_1(j,k)==0
%             for i=1:size(img1,3)
%                 img1_temp(j,k,i)=255;
%             end
            img1_temp(j,k,1)=0;
            img1_temp(j,k,2)=0;
            img1_temp(j,k,3)=255;
        end
    end
end
imwrite(uint8(img1_temp),[img_path,'results_',num2str(para_energy(1)),'_',num2str(para_energy(2)),'/',img_name,'_seg1_post.bmp'],'bmp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accu=0;