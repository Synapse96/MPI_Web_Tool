function [IP,IN,ipos,ineg,maPos,maNeg,envP,envN] = constructEnvelope(I,th_in)
% Construct waveform envelope signal 
%
%
% Description:
%    Threshold is the value that is used to compare with each pixel
%    intensity in each vertical line, the pixel furthest from the basline
%    with intensity higher than the threshold is considered as the envelope
%    in its own vertical line;
%    
%    Note that the baseline is 64 (middle horizontal line in the image) and 
%    it can vary (preferablly allow user to define) 

xaxisL = size(I,2);
w = round(size(I,1)/2);
IP = I(1:w,:);
IN = I((w+1):end,:);

if ~isa(IP,'double')
    ipos = sum(IP);
    ineg = sum(IN);
else
    ipos = std(IP);
    ineg = std(IN);
end

maPos = median(ipos);
maNeg = median(ineg);  

th_orig = th_in;     

for i = 1:size(IP,2)
    rP = IP(:,i);
    rP(rP<th_orig) = 0;
    idP = find(rP>0);
    if isempty(idP)
        envP(i) = 0;
    else
        envP(i) = w-min(idP);
    end
    
    rN = IN(:,i);
    rN(rN<th_orig) = 0;
    idN = find(rN>0);
    if isempty(idN)
        envN(i) = 0;
    else
        envN(i) = max(idN);
        
        if envN(i)>60 && envN(i)<=64
            if sum(rN==0)>= 32
                envN(i)=0;
            end
        end
        
        
    end
end
