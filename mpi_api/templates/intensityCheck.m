function [ tpNew]=intensityCheck(tp,intensity,indicator)

tp1 = tp-5;
tp2 = tp+5;

if tp1 <=0
    tp1 =1;
end

id = [tp1 tp2];
tpNew=zeros(1,length(tp));
tpNew = tpNew';

for i = 1:length(tp)
    
    blk = intensity(id(i,1):id(i,2));
    
    [val(i), k]=max(blk);
    
    tpNew(i)= id(i,1)+k-1;
    
    if nargin <3
        if tpNew(i)-tp(i)>=4
            tpNew(i)=tp(i);
        end
    end
end








