function phi1=get_force_new(img1,idx,upts,vpts,phi,alpha1,alpha2,alpha3,alpha4,alpha5,alpha7,nn)
% img1{nn}=colorspace('HSI<-RGB',img1{nn});
img{1,1}=double(img1{nn}(:,:,1));
img{1,2}=double(img1{nn}(:,:,2));
img{1,3}=double(img1{nn}(:,:,3));
% alpha22=0.05;
idx1=idx{nn};
vpts1=vpts{nn};
%% exterior region
for i=1:size(idx1,1)
    k1=vpts1(find(abs(img{1,1}(idx1(i))-img{1,1}(vpts1))<=alpha2));
    k2=k1(find(abs(img{1,2}(idx1(i))-img{1,2}(k1))<=alpha2));
    k3=k2(find(abs(img{1,3}(idx1(i))-img{1,3}(k2))<=alpha2));
%  k3=vpts1(find(abs(img{1,2}(idx1(i))-img{1,2}(vpts1))<=alpha22));
    E2(i)=size(k3,1)/size(vpts1,1);
end
E2=E2';
%%
kk=1;
phi1=phi{nn};
E=[];
for ki=1:length(img1)
    E1=[];
    %     E2=[];
    if ki==nn
        continue;
    end
%     img1{ki}=colorspace('HSI<-RGB',img1{ki});
    img{2,1}=double(img1{ki}(:,:,1));
    img{2,2}=double(img1{ki}(:,:,2));
    img{2,3}=double(img1{ki}(:,:,3));
    upts2=upts{ki};
    
    for i=1:size(idx1,1)
        k1=upts2(find(abs(img{1,1}(idx1(i))-img{2,1}(upts2))<=alpha2));
        k2=k1(find(abs(img{1,2}(idx1(i))-img{2,2}(k1))<=alpha2));
        k3=k2(find(abs(img{1,3}(idx1(i))-img{2,3}(k2))<=alpha2));
%  k3=upts2(find(abs(img{1,2}(idx1(i))-img{2,2}(upts2))<=alpha22));
        E1(i)=size(k3,1)/size(upts2,1);
    end
    if isempty(upts2)
        E1=zeros(1,length(idx1));
    end
    E=[E,E1'];
    kk=kk+1;
end

%%
kk=1;
phi1=phi{nn};
E3=[];
for ki=1:length(img1)
    E1=[];
    %     E2=[];
    if ki==nn
        continue;
    end
%     img1{ki}=colorspace('HSI<-RGB',img1{ki});
    img{2,1}=double(img1{ki}(:,:,1));
    img{2,2}=double(img1{ki}(:,:,2));
    img{2,3}=double(img1{ki}(:,:,3));
    vpts2=vpts{ki};
    
    for i=1:size(idx1,1)
        k1=vpts2(find(abs(img{1,1}(idx1(i))-img{2,1}(vpts2))<=alpha2));
        k2=k1(find(abs(img{1,2}(idx1(i))-img{2,2}(k1))<=alpha2));
        k3=k2(find(abs(img{1,3}(idx1(i))-img{2,3}(k2))<=alpha2));
%  k3=vpts2(find(abs(img{1,2}(idx1(i))-img{1,2}(vpts2))<=alpha22));
        E1(i)=size(k3,1)/size(vpts2,1);
    end
    if isempty(vpts2)
        E1=zeros(1,length(idx1));
    end
    E3=[E3,E1'];
    kk=kk+1;
end
%% interior force
kk=1;
beta=zeros(length(alpha3)-1,1);
gamma=beta;
for i=1:length(alpha3)
    if i~=nn
        beta(kk)=alpha3(i);
        gamma(kk)=alpha7(i);
        kk=kk+1;
    end
end
if  isempty(vpts1)
    F1=zeros(size(E,1),1);
else
    F= E*beta+alpha4*E2+E3*gamma;
    F=F./max(abs(F));%normalization
    F1=F;
end
%%
curvature = get_curvature(phi1,idx1);  %  curvature penalty
nei=zeros(length(F1),1);
nei=alpha5+nei;
dphidt = F1 + alpha1*curvature + nei;  % gradient descent to minimize energy
%-- maintain the CFL condition
dt = .45/(max(dphidt)+eps);
%-- evolve the curve
phi1(idx1) = phi1(idx1) + dt.*dphidt;
%-- Keep SDF smooth
phi1 = sussman(phi1, .5);



function showCurveAndPhi(I, phi, i)
imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
% contour(phi, [0 0], 'g','LineWidth',4);
contour(phi, [0 0], 'k','LineWidth',2);
hold off; title([num2str(i) ' Iterations']); drawnow;

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