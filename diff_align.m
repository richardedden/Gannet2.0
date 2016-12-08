function [out_v1 out_v2] = diff_align(vector1, vector2, range)

cor_vec=zeros(1,21);
vec1=real(vector1(range));
vec2=(vector2(range));
vec1=[zeros(1,11) vec1 zeros(1,11)];
vect2=[zeros(1,11) vec2 zeros(1,11)];
%figure(99)
%plot(real([vec1.' vect2.']))
%figure(98)
for phase=1:20
    vec2=real(vect2*exp(1i*pi/180*(phase-10)));
for ii=1:21
    %In the future, we might want to edit this to be robust to other
    %resultions of the spectrum....
    shift=ii-11;
    thing=circshift(vec2,[1 shift]);
    diff=vec1-thing;
    buff=20;
    out=std(diff(buff:(end-buff)));
    cor_vec(ii)=out;

end
    
%size(cor_vec)
%figure(20)
%plot(1:21,cor_vec)
[number index]=min(cor_vec,[],2);
phase_vec(phase)=number;
index_vec(phase)=index;
end
[number2 index2]=min(phase_vec,[],2);

out_v1=vector1;
out_v2=circshift(vector2*exp(1i*pi/180*(index2-10)),[1 index_vec(index2)-11]);
end