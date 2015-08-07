function vec1=get_indx_super(phi,label,dd)
vec1=zeros(1,length(dd));
% vec2=vec1;
width=1.2;
mask1=find(phi<-width);
% mask2=find(phi>width);
mask3=find(phi <= width & phi >= -width);

vec1(label(mask1))=1;
vec1(label(mask3))=-1;