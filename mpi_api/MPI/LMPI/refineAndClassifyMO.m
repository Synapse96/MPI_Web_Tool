function [tp2newv2,MOindict]=refineAndClassifyMO(tempIN,turningPoints,envelope,N,s6v1)
% This function refines detected MO using intensity

tempIN2 = tempIN./envelope;
tempIN2(tempIN2==inf)=0;tempIN2(tempIN2==-inf)=0;
tp2newv2 = turningPoints;
k = 10;
MOindict = nan(1,length(turningPoints));

for i = 1:length(turningPoints)
    p = turningPoints(i); fprintf('MO: %d ',p);
    tempid = turningPoints(i)-k:turningPoints(i)+k;
    
    tempsigGx = tempIN2(tempid);
    [peak,locs] = max(tempsigGx);
    locsOrig = locs;
    thPeaks = 0.05;
    
    [peak3,locs3] = max(tempIN(tempid));
    locs3 = tempid(locs3);
    
    pnew = tempid(locs);
    if pnew+k<=N
        tempid2 = pnew:pnew+k;
    else
        tempid2 = pnew:N;
    end
    
    [minval,minid] = min(tempIN2(tempid2));
    minid = tempid2(minid);
    [peak2,locs2] = findpeaks(tempsigGx);
    locs2 = tempid(locs2);
    % remove non significant local maxima 
    locs2(peak2<0.15|envelope(locs2)<5|s6v1(locs2)<0.5|s6v1(locs2)<median(s6v1)/5)=[];
    locs2 = locs2(locs2<pnew);
    peak2 = tempIN2(locs2);
    fprintf('[%d %d] (%0.2f %0.2f) ',pnew,minid,peak ,minval);
    if locs3-pnew==1
        pnew = locs3;
    end
    
    if peak<thPeaks
        fprintf('Very small value, might discard %d ',p);
        tp2newv2(i)=0;MOindict(i)=inf;
    else
        if pnew>p && length(locs2)==1 && locs2<p
            tp2newv2(i)=locs2+round(diff([locs2 p])/2);
            fprintf(2,'relatively good ');
            MOindict(i)=1;
        elseif isempty(locs2)&& peak>=1
            MOindict(i)=0; fprintf(2,'good (might be slightly thick or tilt or merge with flow)');
            if  minid-pnew<=6
                tp2newv2(i)= pnew+1;
            elseif minid-pnew>6 && minid-pnew<9
                tp2newv2(i)= pnew+2;
            elseif minid-pnew>=9
                if  tempIN(minid)<-20&& tempIN(pnew)>20
                    tp2newv2(i)= pnew+3;
                else
                    tp2newv2(i)= pnew+2;
                end
            end
        elseif length(locs2)==1&&peak>=1
            if peak2>peak*0.5 || p-pnew<4
                tp2newv2(i)= locs2+round(diff([locs2 pnew])/2);MOindict(i)=1;
            else
                tp2newv2(i)= pnew+1;MOindict(i)=0;
            end
        elseif length(locs2)>1 && peak>=1
            if (peak2(end)>peak*0.5 || p-pnew<4)&&peak>=1.5
                tp2newv2(i)= locs2(end)+round(diff([locs2(end) pnew])/2);MOindict(i)=1;
            else
                tp2newv2(i)= pnew+1;MOindict(i)=0;
            end
        else
            fprintf(2,'Not very distinct or bright; No. locs2 %d; pnew %d(%1.2f) and minid %d(%1.2f) ',length(locs2),pnew,peak,minid,minval);
            if length(locs2)==1 && pnew-locs2<=5
                pnew = locs2+1;
            end
            if peak>0.5&&minval<0&&p>pnew
                tp2newv2(i) = pnew+1;MOindict(i)=2;
            elseif peak>0.5&&minval>0
                tp2newv2(i) = pnew+2;MOindict(i)=2;
            elseif pnew == p
                tp2newv2(i) = p;MOindict(i)=2;
            elseif peak>0.4&&minval<0&&p>pnew
                tp2newv2(i) = pnew+1;MOindict(i)=2;
            else
                fprintf('MO is not Not good  ');MOindict(i)=5;
            end
        end
    end
    if tp2newv2(i)>p+5
        tp2newv2(i)=p;
    end
    fprintf('\n')
end
MOindict = MOindict(tp2newv2>0);
tp2newv2 = tp2newv2(tp2newv2>0);

