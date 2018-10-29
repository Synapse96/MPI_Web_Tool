function [ao,indict] = refineAndClassifyAO(tempid,sumGx1,sumGx1Orig,N,p,s1,envP1)
% determine the quality of detected AO by check the intensity variation as
% well as shape (envelope ratio) 

ao=0;
indict = 0;
w = 5;
k = 10;
tempsigGx = sumGx1(tempid);
[peak,locs] = max(tempsigGx);
thPeaks = 0.05;

[peak2,locs2] = findpeaks(tempsigGx);
locs2(peak2<0.1)=[];
peak2(peak2<0.1)=[];
locs2 = tempid(locs2);

if p+10<=N && s1(p+1)/max(s1(p:p+1+10))<0.75
    peak2(p-locs2>=(w+k)-1)=[];
    locs2(p-locs2>=(w+k)-1)=[];    
end

pnew = tempid(locs);
if pnew+k<=N
    tempid2 = pnew:pnew+k;
else
    tempid2 = pnew:N;
end

[minval,minid] = min(sumGx1(tempid2));
minid =tempid2(minid);

[peak3,locs3] = max(sumGx1Orig(tempid));
locs3 = tempid(locs3);
if locs3-pnew==1 && peak<1.5
    pnew = locs3;
end

if peak<thPeaks
    fprintf('Very small value, might discard %d\n',p);
    ao=0;
    indict =inf;
end

if (peak>=2 && minval<-0.5&& (minid-pnew<=w+1 &&minid>pnew))
    fprintf('C1: good aortic opening %d [%d %d]\n',p,pnew,minid);
    ao= pnew+2;
elseif (peak>=2 && minval<-0.5&& (minid-pnew<k &&minid>pnew+w+1))
    fprintf('C1: good thick aortic opening %d [%d %d]\n',p,pnew,minid);
    ao= pnew+3;
elseif peak>=1 && minval<-0.25&& (minid-pnew<=w+1 &&minid>pnew)
    fprintf('C1: good aortic opening (not perfect) %d [%d %d]\n',p,pnew,minid);
    ao= pnew+1;
elseif peak<1
    fprintf('C2: suboptimal aortic opening %d [%d %d]\n',p,pnew,minid);
    if minval>0
        indict = 3;
    else
        indict = 1;
    end
    locs2 = locs2(locs2<pnew);
    peak2 = peak2(locs2<pnew);
    if length(locs2)==1 && locs2+w<=N && locs2-w>0
        if mean(s1(locs2-5:locs2-1))>s1(locs2)&&mean(s1(locs2+1:locs2+5))>s1(locs2)
            locs2=[];
            peak2=[];
        end
    end
    
    if length(locs2)==1 && peak2>=1/2*sumGx1(pnew)
        ao = locs2+1;
    else
        indict = indict+1;
        if pnew+1+10<=N
            envRatio = s1(pnew+1)/max(s1(pnew:pnew+k+1));
            envRatio2 = envP1(pnew+1)/max(envP1(pnew:pnew+k+1));
        else envRatio = 0.5;
        end
        if envRatio <0.7
            ao= pnew+1;
        elseif envRatio<1 && envRatio>=0.7
            if envRatio<0.8 && envRatio2<0.7
                ao= pnew+1;
            elseif pnew-10>0 && (min(sumGx1(pnew-k:pnew))<-2||min(sumGx1(pnew-k:pnew))==0)
                ao= pnew+1;
            else
                ao=0; indict =indict+inf;
                fprintf('Not good , envRatio %d\n',envRatio);
            end
            
        end
    end
    
elseif peak<0.2 && minval<0.1
    fprintf('C3: Cannot find a good opening %d [%d %d]\n',p,pnew,minid);
    ao=0;indict = inf;
    
else fprintf('Check this aortic opening %d\n',p);
    if minval>0
        indict = 3;
    else
        indict = 1;
    end
    
    if (peak>=3 && minval<-0.25&& (minid-pnew==10))
        fprintf('C1v2: good thick aortic opening %d [%d %d]\n',p,pnew,minid);
        ao= pnew+3;indict = 2;
    elseif peak>=1&& pnew+k+1<=N
        envRatio = s1(pnew+1)/max(s1(pnew:pnew+k+1));
        locs2 = locs2(locs2<pnew);
        peak2 = peak2(locs2<pnew);
        
        if sum(pnew-locs2>w)~=0
            if sum(sumGx1(locs2(pnew-locs2>w):pnew)<-1)~=0
                peak2(pnew-locs2>w) = []; locs2(pnew-locs2>w) = [];
            end
        end
        
        if length(locs2)==1 && peak2>1/2*peak && pnew-locs2>w
            fprintf(2,'check here for %d\n',pnew);
            pnew = locs2-round(diff([pnew locs2])/2);
            indict =indict+1;
        end
        if envRatio > 0.3&& envRatio<0.7
            if p==pnew
                ao=pnew+1;
            else
                if minval>0
                    ao= pnew+3;
                elseif minid-pnew>=6
                    ao= pnew+2;
                else
                    ao= pnew+1;
                end
            end
            indict =indict+1;
            fprintf('C2v2: suboptimal aortic opening %d [%d %d]\n',p,pnew,minid);
        elseif envRatio<0.8
            if minval>0
                ao= pnew+2;
            else
                ao= pnew+1;
            end
            indict =indict+1;
        end
    else
        fprintf('NOT SURE: max min maxid minid\n [%d %d %d %d]\n', peak,minval,pnew,minid)
        indict =indict+inf;
    end
end

