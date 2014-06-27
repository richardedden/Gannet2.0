
function Output = image_center(Input,Window_size)
    %
    Output=zeros([Window_size Window_size]);
    size1=size(Input,1)-1;
    size2=size(Input,2)-1;
    offset1=round((Window_size-size1)/2);
    offset2=round((Window_size-size2)/2);
    Output(offset1:(offset1+size1),offset2:(offset2+size2))=Input;
end
