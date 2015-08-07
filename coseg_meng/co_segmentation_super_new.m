
function phi = co_segmentation_super_new(img,mm,n_for_show,num_iter,super,is_show,parameter,large_setp)
%%
% parameters;
width=parameter.width;
alpha1=parameter.alpha1;
alpha2=parameter.alpha2;
alpha3=parameter.alpha3;
alpha4=parameter.alpha4;
alpha5=parameter.alpha5;
alpha7=parameter.alpha7;


%-- initialization
phi=[];
for i=1:length(mm)
    phi{i} = mask2phi(mm{i});
    super(i).vec=get_indx_super(phi{i},super(i).label,super(i).num);
    % phi2 = mask2phi(m2);
end

%% iteration
flag=ones(1,length(phi));% 1: image with evolution. 0: image without evolution
ratio=ones(1,length(phi));
for its = 1:num_iter % number of iteration
    
    disp([num2str(its),'...',num2str(num_iter)]);
    idx=[];
    upts=[];
    vpts=[];
    for i=1:length(phi)
        idx{i}=find(phi{i} <= width & phi{i} >= -width);
        upts{i} = find(phi{i}<=0);
        vpts{i} = find(phi{i}>0);
    end
    %% evolution
    for i=1:length(phi)
        if isempty(idx{i})% if there are no foreground region
            flag(i)=0;
        end
        if flag(i)
            if its==1
                phi_ori=phi;
            end
            if (mod(its,10) == 1) && its~=1
                [phi,flag(i),ratio(i)]=rescale_phi(phi,phi_ori,super,i); % refine the curve based on the superpixel per 10 iterations
                phi_ori{i}=phi{i};
            end
            if is_show==1
                if(mod(its,n_for_show) == 1)
                    figure(i);
                    showCurveAndPhi(uint8(img{i}),phi{i},its);
                    % % %             saveas(figure(i),['.\res\super_',num2str(i),'_',num2str(its),'_',num2str(mengi),'_image'],'jpg');
                end
            end
            %% run the co-segmentation on the its-th image
            phi{i}=get_force_super_new(img,idx,phi,alpha1,alpha2,alpha3,alpha4,alpha5,alpha7,i,super,large_setp);
        end
    end
    %% foreground superpixel based on new segment
    for i=1:length(mm)
        super(i).vec=get_indx_super(phi{i},super(i).label,super(i).num);
    end
    
end


%-- Displays the image with curve superimposed
% % % function showCurveAndPhi(I, phi, i)
% % % imshow(I,'initialmagnification','fit','displayrange',[0 255],'Border','tight'); 
% % % set (gcf,'Position',[300,300,500,500*(size(I,1)/size(I,2))]);
% % % axis normal;
% % % hold on;
% % % contour(phi, [0 0], 'g','LineWidth',20);
% % % contour(phi, [0 0], 'k','LineWidth',10);
% % % hold off; %title([num2str(i) ' Iterations']); 
% % % drawnow;


function showCurveAndPhi(I, phi, i)
imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
contour(phi, [0 0], 'g','LineWidth',4);
contour(phi, [0 0], 'k','LineWidth',2);
hold off; 
title(['Image ',num2str(i)]);
% title(['Image ',num2str(i) 'of Iteration ',num2str(j)]); 
drawnow;

%-- converts a mask to a SDF
function phi = mask2phi(init_a)
phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;

%-- compute curvature along SDF
function curvature = get_curvature(phi,idx)
[dimy, dimx] = size(phi);
[y x] = ind2sub([dimy,dimx],idx);  % get subscripts

%-- get subscripts of neighbors
ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

%-- bounds checking
ym1(ym1<1) = 1; xm1(xm1<1) = 1;
yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;

%-- get indexes for 8 neighbors
idup = sub2ind(size(phi),yp1,x);
iddn = sub2ind(size(phi),ym1,x);
idlt = sub2ind(size(phi),y,xm1);
idrt = sub2ind(size(phi),y,xp1);
idul = sub2ind(size(phi),yp1,xm1);
idur = sub2ind(size(phi),yp1,xp1);
iddl = sub2ind(size(phi),ym1,xm1);
iddr = sub2ind(size(phi),ym1,xp1);

%-- get central derivatives of SDF at x,y
phi_x  = -phi(idlt)+phi(idrt);
phi_y  = -phi(iddn)+phi(idup);
phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
    +0.25*phi(iddr)+0.25*phi(idul);
phi_x2 = phi_x.^2;
phi_y2 = phi_y.^2;

%-- compute curvature (Kappa)
curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
    (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);

%-- Converts image to one channel (grayscale) double
function img = im2graydouble(img)
[dimy, dimx, c] = size(img);
if(isfloat(img)) % image is a double
    if(c==3)
        img = rgb2gray(uint8(img));
    end
else           % image is a int
    if(c==3)
        img = rgb2gray(img);
    end
    img = double(img);
end

%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
% forward/backward differences
a = D - shiftR(D); % backward
b = shiftL(D) - D; % forward
c = D - shiftD(D); % backward
d = shiftU(D) - D; % forward

a_p = a;  a_n = a; % a+ and a-
b_p = b;  b_n = b;
c_p = c;  c_n = c;
d_p = d;  d_n = d;

a_p(a < 0) = 0;
a_n(a > 0) = 0;
b_p(b < 0) = 0;
b_n(b > 0) = 0;
c_p(c < 0) = 0;
c_n(c > 0) = 0;
d_p(d < 0) = 0;
d_n(d > 0) = 0;

dD = zeros(size(D));
D_neg_ind = find(D < 0);
D_pos_ind = find(D > 0);
dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
    + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
    + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;

D = D - dt .* sussman_sign(D) .* dD;

%-- whole matrix derivatives
function shift = shiftD(M)
shift = shiftR(M')';

function shift = shiftL(M)
shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
shift = shiftL(M')';

function S = sussman_sign(D)
S = D ./ sqrt(D.^2 + 1);






