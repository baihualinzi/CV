function mm=run_coseg(img_path, para_energy, is_show,ratio,iter)

mkdir([img_path,'results']);% generate the folder for output
img_name=dir([img_path,'*.jpg']);% the name of the images

%% pre-processing of superpixel generation
%%
tic
img=[]; % save image
super=[]; % save superpixel
out='.\image\SLIC\'; % for SLIC code
mkdir(out); 
func='SLICSuperpixelSegmentation';
cand_param=[10,1000]; % 1000 is the number of superpixel
%% pre-processing
disp('Pre-processing of SLIC over-segmentation...');
for i=1:length(img_name)
    disp([num2str(i),'...',num2str(length(img_name))]);
    img{i}=double(imread([img_path,img_name(i).name]));% read the image
    %%
    imwrite(uint8(img{i}),[img_path,strrep(img_name(i).name,'.jpg','.bmp')]);% for SLIC
    file=[img_path,strrep(img_name(i).name,'.jpg','.bmp')];
    commond=[func,' ',file,' ',num2str(cand_param(1)),' ',num2str(cand_param(2)),' ',out];
    system(commond);% SLIC
    label_idx=fopen([out,strrep(img_name(i).name,'.jpg','.dat')],'r');
    label=(fread(label_idx,[size(img{i},2),size(img{i},1)],'int'))';
    fclose(label_idx);
    label=label+1; % the superpixels of SLIC
    dd=[];
    color=[];
    % save each superpixel into super
    for j=1:max(label(:))
        ind=find(label==j);
        a1=img{i}(:,:,1);
        a2=img{i}(:,:,2);
        a3=img{i}(:,:,3);
        R=mean(a1(ind));
        G=mean(a2(ind));
        B=mean(a3(ind));
        dd(j)=length(find(label==j));
        color(j,:)=[R,G,B];
    end
    super(i).label=label;
    super(i).num=dd;
    super(i).color=color;
end

%% initial curves mm
mm=cell(1,length(img_name));
ra=0.12;% the width of the distance between the curve and the image edge.

for  j=1:length(img_name)
%     delt=40;
%  mm{j}(delt+1:end-delt,delt+1:end-delt)=1;%  mm{j}=zeros(size(img{j},1),size(img{j},2));
%  mm{j}(delt:size(img{j},1)-delt,delt:size(img{j},2)-delt)=1

 delt=round(min(size(img{j},1),size(img{j},2))*ra);
% delt=80;
  mm{j}=zeros(size(img{j},1),size(img{j},2));
  mm{j}(delt+1:end-delt,delt+1:end-delt)=1;
end


%%set parameters in the co-segmentation
parameter=[];
parameter.width=1.5;%for level-set band determinition
parameter.alpha1=0.001; %curvature
parameter.alpha2=18; %color band
parameter.alpha3=-para_energy(1)*(ones(length(img_name),1)/(length(img_name)-1)); % for inner region
parameter.alpha4=5; % for exterior region
parameter.alpha5=0.001;%for \nu
parameter.alpha7=para_energy(2)*(ones(length(img_name),1)/(length(img_name)-1)); % background term

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial co-segmentation
mm_temp=[];
img_temp=[];
for j=1:length(img)
    img_temp{j} = imresize(img{j},ratio(1),'nearest');  %-- resize the image into small size
    mm_temp{j} = imresize(mm{j},ratio(1),'nearest');  %     for fast computation
end

phi = co_segmentation_new(img_temp,mm_temp,10,iter(1),parameter); %-- Run co-segmentation
%     ratio_1=get_save_results(img_path,phi1,img1,m1_t);
mm=[];
for j=1:length(phi)
    mm_s=zeros(size(phi{j},1),size(phi{j},2));
    mm_s(find(phi{j}<=0))=1;
    mm_s=post_processing(mm_s);
    mm_s=imresize(mm_s,[size(img{j},1),size(img{j},2)],'nearest');
    mm{j}=mm_s;
end

%% super based cosegmentation
for i=2:size(ratio,2)
    if i==1
        large_setp=1;
    else
        large_setp=0;
    end
%     save meng_12q meng_12q;
    %%scale the images and the maskes
    mm_temp=[];
    img_temp=[];
    for j=1:length(img)
        img_temp{j} = imresize(img{j},ratio(i),'nearest');  %-- resize the image into small size
        mm_temp{j} = imresize(mm{j},ratio(i),'nearest');  %     for fast computation
    end
    %%scale the superpixel information
    super_temp=[];
    for j=1:length(super)
        super_temp(j).label=imresize(super(j).label,ratio(i),'nearest');
        super_temp(j).num=round(super(j).num*ratio(i));
        super_temp(j).color=super(j).color;
    end
%     mengi=i;
%     save mengi mengi;
    disp(['The ',num2str(i),'-th hierachy of ', num2str(length(ratio)),' (',num2str(ratio(i)),')']);
    phi = co_segmentation_super_new(img_temp,mm_temp,10,iter(i),super_temp,is_show,parameter,large_setp); %-- Run co-segmentation
    %% new segment
    mm=[];
    for j=1:length(phi)
        mm_s=zeros(size(phi{j},1),size(phi{j},2));
        mm_s(find(phi{j}<=0))=1;
        mm_s=post_processing(mm_s);
        mm_s=imresize(mm_s,[size(img{j},1),size(img{j},2)],'nearest');
        mm{j}=mm_s;
    end
end
%%
mkdir([img_path,'results_',num2str(para_energy(1)),'_',num2str(para_energy(2))]);
for i=1:length(phi)
    get_save_results(img_path,phi{i},img_temp{i},para_energy,strrep(img_name(i).name,'.jpg',''));
end