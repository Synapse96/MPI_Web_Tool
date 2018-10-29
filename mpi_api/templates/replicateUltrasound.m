function Im =  replicateUltrasound (I,mark)
% for display purpose 

I = imadjust(I,[0 1],[0 0.8],0.6);

th = graythresh(I);

I = imadjust(I,[0.05 1],[0 0.8],0.5);

if nargin ==1
    [numr, numc] = size(I);
    I = imresize(I, [floor(numc/2.7) numc]);
    Im = I;
    
else
    mask = I<0;
    mask(:,mark) = 1;
    addmark = imoverlay(I,mask,[0 1 0]);
    Im = imresize(addmark,[210 560]);
    
end



end