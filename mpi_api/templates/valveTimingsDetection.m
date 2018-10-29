function [AO,valueAO,AC,valueAC,MO,valueMO,MC,valueMC,s5v2,sumtempIP,tempIN,MCorig,sumGx1,sumGx1Orig]=valveTimingsDetection(data)
% Find valve events 
%
%
% Description:
%    This is the initial attempt of detecting valve events using envelope
%    signals, and the result of which might contain wrongly detected events
%    or miss some events


I_orig = data;

I = imadjust(data);
I2 = I;
I = medfilt2(I,[5 5]);

[counts,~]=imhist(I);
[level1]=triangle_th(counts,256);

[counts2,~]=imhist(I_orig);
[level2]=triangle_th(counts2,256);

[IP,IN,ipos,ineg,maPos,maNeg,envP,envN] = constructEnvelope(I,level1);
[IP1,IN1,ipos1,ineg1,maPos1,maNeg1,envP1,envN1] = constructEnvelope(I_orig,level2);
[IP2,IN2,ipos2,ineg2,maPos2,maNeg2,envP2,envN2] = constructEnvelope(I2,level2);

[Gx1, Gy1] = imgradientxy(imadjust(data));
sumGx1 = sum(Gx1(1:64,:));
sumGx1Orig = sumGx1;
tempImg = Gx1.*imadjust(data);
tempIP = tempImg(1:64,:);
sumtempIP = sum(tempIP);
sumtempIN = sum(tempImg(65:128,:));

envN = movavgFilt(envN,5,'Center');
envP = movavgFilt(envP,5,'Center');
s5 = movavgFilt(sum(IP),5,'Center');
s6 = movavgFilt(sum(IN),5,'Center');


s1 = envP;
s2 = envP.*ipos;
s3 = envN;
s4 = envN.*ineg;

s5v2 = sum(IP2); 
s6v1 = sum(IN1);
N = size(I,2);
k1 = 10;
k2 = 20;
w = 5;
k=k1;


sumGx1 = sumGx1./s1;
sumGx1(sumGx1==inf)=0;sumGx1(sumGx1==-inf)=0;


%% Mitral closure
[MC,comp4V2,MCorig] = findMC(N,k1,k2,s1,s3,s4,w,IP1,IP2,IN2,envP,s5v2,sumtempIP,s5);

%% Aortic opening
[compAOraw,~,compAO]= getCompSigForAO(N,k2,w,k1,s1,s3,s5,s6);
cprintf('Hyperlinks','ValveTiming function: AO\n');

if sum(compAO)==0
    AO = [];
    fprintf('Composite signal is zero, upper wave not obvious\n');
else
    
    thCompAO = median(compAO(compAO~=0))/2;
    tpCompAO = turning_points(compAO,thCompAO);
    tpCompAOpos = tpCompAO(tpCompAO>0);
    tpCompAOposOrig = tpCompAOpos;
    for i = 1:length(tpCompAOpos)       
        tpCompAOpos(i)=intensityCheck(tpCompAOpos(i),sum(IP1));        
    end
    tpCompAOpos;
    
    clear i
    for i = 1:length(tpCompAOpos)-1
        p1 = tpCompAOpos(i);
        p2 = tpCompAOpos(i+1);
        
        if p2-p1<100 && p2-p1>=k1*3
            if  compAOraw(p1)>compAOraw(p2) 
                tpCompAOpos(i+1)=0;
            else
                tpCompAOpos(i)=0;
            end
            
        elseif p2-p1<k1*3
             if min(s1(p1:p2))==s1(p1)
                if sum(compAO(p1:p2)==0)~=0
                    tpCompAOpos(i+1)=0;
                elseif abs(compAO(p1)-compAO(p2))/min([compAO(p1) compAO(p2)])<0.5*min([compAO(p1) compAO(p2)]) && p2-p1<k1
                    tpCompAOpos(i+1)=0;
                elseif p2-p1<k1
                    tpCompAOpos(i)=0;
                elseif min([s1(p1)/max(s1(p1:p1+k1+w)) s1(p2)/max(s1(p2:p2+15))])== s1(p1)/max(s1(p1:p1+k1+w))&&...
                        s1(p1)/max(s1(p1:p1+k1+w))<0.6 && s1(p1)/max(s1(p1:p1+k1+w))>0.2
                    tpCompAOpos(i+1)=0;
                elseif min([s1(p1)/max(s1(p1:p1+k1+w)) s1(p2)/max(s1(p2:p2+k1+w))])== s1(p1)/max(s1(p1:p1+k1+w))&&...
                        s1(p1)/max(s1(p1:p1+k1+w))<0.7 && s1(p1)/max(s1(p1:p1+k1+w))>0.6
                    tpCompAOpos(i+1)=0;
                end
             else
                 tpCompAOpos(i)=0;
            end        
        end        
    end
    

    tpCompAOpos=tpCompAOpos(tpCompAOpos>0);
    tpCompAOposRefined=tpCompAOpos;
    
    clear i p p1 p2
    
    tpCompAOFinal = tpCompAOposRefined;
    AOindict = zeros(1,length(tpCompAOFinal));
    plotid = 1;

    for i = 1:length(tpCompAOposRefined)
        p = tpCompAOposRefined(i);
        
        pOrig = tpCompAOposOrig(abs((p-tpCompAOposOrig))==min(abs((p-tpCompAOposOrig))));
        p1 = tpCompAOposRefined(i)-k2;
        tempid = p1:pOrig;
        p1v2=find(compAO(tempid)==0); 
        
        tempidend = p:tpCompAOposRefined(i)+k1;
        tempend=find(compAO(tempidend)==0);
        if ~isempty(tempend)
            tempend=tempidend(tempend);tempend = tempend(1);
        else tempend = p+k1;
        end
        
        
        if ~isempty(p1v2)
            p1v2 = p1v2(end);
            tempid =p1+p1v2-1:tempend ;
        else tempid = p1:tempend;
        end
        
        if ~isempty(tempid) && length(tempid)<k1
            tempsig = sum(IP1);
            tempsig = tempsig(tempid);
            if max(tempsig)== tempsig(p-tempid(1)+1) && p+k2<=N
                if s1(p)/max(s1(p:p+k2))<0.3
                    tempid = tempid(1)+w:p+w;
                end
            end
        end
             
        if length(tempid)==3
            tempid = tempid(1)-2:tempid(end); 
            fprintf('\nExtend tempid for %d ', p);
        elseif length(tempid)<3 && ~isempty(tempid)
            
            if sumGx1(tempid(1))<sumGx1(tempid(1)+1)
                tempid =  tempid(1):tempid(end)+(w-1);
            elseif sumGx1(tempid(1))<sumGx1(tempid(1)-1)
                tempid =  tempid(1)-(w-1):tempid(end);
            elseif tempid(1)-k2>0
                [minv,minidx]=min(sumGx1(tempid(1)-k2:tempid(end)));
                if minv<-2
                    tempid = tempid(1)-k2+minidx+2:tempid(end);
                end
            end
        end
        
 
        if length(tempid)>=3
               
            img = data(1:64,tempid);
            [Gmag,Gdir,Gx,Gy]= plotGmag(img,tpCompAOposRefined,tempid,plotid,'no plot',s5v2);
            tempsig2 = sum(Gmag);

            plotid = plotid+1;
            
            if ((max(tempsig2)-min(tempsig2))/min(tempsig2))<1 && ...
                    (s1(tpCompAOFinal(i))/max(s1(tpCompAOFinal(i)+1:tpCompAOFinal(i)+k1+w))>0.6 &&...
                    s1(tpCompAOFinal(i))/max(s1(tpCompAOFinal(i)+1:tpCompAOFinal(i)+k1+w))<1)
                
                fprintf('\n%d std of sum(Gmag) is %f, less than 1,\nnot a good rough area, disgard?',...
                    p,(max(tempsig2)-min(tempsig2))/min(tempsig2));
                tpCompAOFinal(i) = 0; 
                continue
                
            end
            
            [tpCompAOFinal(i),AOindict(i)]= refineAndClassifyAO(tempid,sumGx1,sumGx1Orig,N,p,s1,envP1);
            
        end
    end
    
    AOindict(tpCompAOFinal<=0)=[];
    tpCompAOFinal = tpCompAOFinal(tpCompAOFinal>0);
    AO = [tpCompAOFinal AOindict'];
end


%% Aortic Closure
compACv1 = zeros(1,N);
compACv2 = zeros(1,N);
compACv3 = zeros(1,N);

for i = (k1+w+1):(N-w-k1)
    
    id1 = (i-w-k1):(i-w-1);
    id2 = (i+w+1):(i+w+k1);
    
    compACv1(i) = s3(i)-mean(s3(id1))-mean(s3(id2))-mean(s1(id2));
    compACv2(i) = s3(i)+2*mean(s1(id1))-4*mean(s1(id2))-2*mean(s3(id1));
    compACv3(i) = s4(i)-2*mean(s4(id1))-mean(s4(id2));
end


compACv1(compACv1<0)=0;
compACv2(compACv2<0)=0;
compACv3(compACv3<0)=0;

compACv1 = compACv1.*compACv2/max(compACv2.*compACv1);

compACv4 = s4;
compACv4(compACv1==0)=0;
compACv4(compACv3==0)=0;

if max(compACv4)>w&&max(compACv4)<k1
    thCompAC = median(compACv4(compACv4~=0))/3;
else
thCompAC = median(compACv4(compACv4~=0))/2;
end
tpCompAC = turning_points(compACv4,thCompAC);
tpCompAC = tpCompAC(tpCompAC>0);

clear i
for i = 1:length(tpCompAC)-1
    
    if tpCompAC(i)-k1-w<=0
        tpCompAC(i)=0;
    else
        if tpCompAC(i+1)-tpCompAC(i)>k1*5 && tpCompAC(i+1)-tpCompAC(i)<150
            if sum(s3(tpCompAC(i)-k1-w:tpCompAC(i)-w))>sum(s3(tpCompAC(i+1)-k1-w:tpCompAC(i+1)-w))
                tpCompAC(i)=0;
            else tpCompAC(i+1)=0;
            end
            
        elseif tpCompAC(i+1)-tpCompAC(i)<=k1*5 && tpCompAC(i+1)-tpCompAC(i)>k1*3
            if sum(s1(tpCompAC(i)-k1-w:tpCompAC(i)-w))>sum(s1(tpCompAC(i+1)-k1-w:tpCompAC(i+1)-w))
                tpCompAC(i)=0;
            else tpCompAC(i+1)=0;
            end
        elseif tpCompAC(i+1)-tpCompAC(i)<=k1*3
            if sum(s1(tpCompAC(i)-k1-w:tpCompAC(i)-w))>sum(s1(tpCompAC(i+1)-k1-w:tpCompAC(i+1)-w))||...
                    sum(s3(tpCompAC(i)-k1-w:tpCompAC(i)-w))<sum(s3(tpCompAC(i+1)-k1-w:tpCompAC(i+1)-w))
                tpCompAC(i+1)=0;
            end
            
        end
    end
    
end

tpCompAC = tpCompAC(tpCompAC>0);
AC = tpCompAC;
AC = intensityCheck(AC,sum(IN));

%% Mitral Opening
compMOv1 = zeros(1,N);
compMOv2 = zeros(1,N);

for i = (k+w+1):(N-w-k)
    
    id1 = (i-k):(i-1);
    id2 = (i+1):(i+k);
    id3 = (i-w-k):(i-w-1);
    id4 = (i+w+1):(i+w+k);
    
    compMOv1(i) = s3(i)+ mean(s3(id2))+mean(s3(id4))-4*mean(s3(id1))-mean(s3(id3))-2*mean(s1(id1))-2*mean(s1(id4));
    compMOv2(i) = ineg(i)-mean(ineg(id3));
    
end
compMOv2(compMOv2<0)=0;
compMOv1(compMOv1<0) = 0;
compMOv1 = compMOv1.*compMOv2;


compMOv3 = s6;
compMOv3(compMOv1==0)=0;
compMOv3 = medfilt1(compMOv3,5);


thCompMO = median(compMOv3(compMOv3~=0))/2; 
tpCompMO = turning_points(compMOv3,thCompMO);
tpCompMO = tpCompMO(tpCompMO>0);

tpCompMO(s1(tpCompMO)==0&s3(tpCompMO)==0)=[];
for j = 1:(length(tpCompMO)-1)
    
    if tpCompMO(j+1)-tpCompMO(j)<=k1*3
        tpCompMO(j)=0;
    elseif tpCompMO(j+1)-tpCompMO(j)>k1*3 && tpCompMO(j+1)-tpCompMO(j)<100
        
        if tpCompMO(j+1)+5<=N && mean(s1(tpCompMO(j)+1:tpCompMO(j)+5))>mean(s1(tpCompMO(j+1)+1:tpCompMO(j+1)+5))&&(mean(s3(tpCompMO(j)+1:tpCompMO(j)+5))<mean(s3(tpCompMO(j+1)+1:tpCompMO(j+1)+5))||...
                (tpCompMO(j+1)+10<=N && mean(s3(tpCompMO(j)+5:tpCompMO(j)+10))<mean(s3(tpCompMO(j+1)+5:tpCompMO(j+1)+10))))
            tpCompMO(j)=0;
        elseif tpCompMO(j+1)+5<=N && mean(s3(tpCompMO(j)+1:tpCompMO(j)+5))<mean(s3(tpCompMO(j+1)+1:tpCompMO(j+1)+5)) && s3(tpCompMO(j))<s3(tpCompMO(j+1)) 
            tpCompMO(j+1)=0;
        end
         
    elseif tpCompMO(j+1)-tpCompMO(j)>100 && tpCompMO(j+1)-tpCompMO(j)<120 && tpCompMO(j+1)+10<=N
        if sum(s3(tpCompMO(j+1):tpCompMO(j+1)+10))>sum(s3(tpCompMO(j):tpCompMO(j)+10))
            tpCompMO(j)=0;
        else
            tpCompMO(j+1)=0;
        end
        
    end
    
    
    
end
tpCompMO = tpCompMO(tpCompMO>0);

tempIN = sum(Gx1(65:128,:));
tempIN2 = tempIN./s3;
tempIN2(tempIN2==inf)=0;
tempIN2(tempIN2==-inf)=0;

cprintf('Hyperlinks','ValveTiming function: MO\n');

[tpCompMO,MOindict]=refineAndClassifyMO(tempIN,tpCompMO,s3,N,s6v1);
MO = [tpCompMO MOindict'];


valueAO = envP(AO(:,1));
valueAC = envN(AC);
valueMO = envN(MO(:,1));
valueMC = envP(MC);





function [MC,compMC,MCorig] = findMC(N,k1,k2,s1,s3,s4,w,IP1,IP2,IN2,envP,s5v2,sumtempIP,s5)

compMCtemp = zeros(1,N);
compMCtempMod = zeros(1,N);
s5v1=sum(IP1);

for i = (k2+w+1):(N-w-2*k2)
    
    id1 = (i-w-k2):(i-w-k1);
    id2 = (i+w+k2):(i+w+2*k2);
    
    id3 = (i-w-k1):(i-w-1);
    id4 = (i+w+1):(i+w+k1);
    
    compMCtemp(i) = s1(i)+4*mean(s1(id2))-2*mean(s3(id4))-2*mean(s1(id3))-s3(i)-2*mean(s3(id3))-mean(s3(id2));
end
compMCtemp(compMCtemp<0)=0;

envDiff = s1-s3;
envDiff(envDiff<0)=0;
compMCtemp([1:k2+w N-w-2*k2+1:end])=envDiff([1:k2+w N-w-2*k2+1:end]);

for i = (k2+w+1):(N-w-k2)  
    
    id1 = (i-w-k2):(i-w-k1);
    id2 = (i+w+k1):(i+w+k2);
    
    id3 = (i-w-k1):(i-w-1);
    id4 = (i+w+1):(i+w+k1);
    
    compMCtempMod(i) = compMCtemp(i)-6*max(compMCtemp(id1));
end

compMCtempMod(compMCtempMod<0)=0;


compMC = sum(IP1);
compMC(compMCtempMod==0)=0;

thCompMC = median(compMC(compMC~=0))/2;
tpCompMC = turning_points(compMC,thCompMC);
tpCompMCpos = tpCompMC(tpCompMC>0);

% 50ms time interval / sampling frequency
intvl = 0.05/0.002; 

% remove extra detections
clear i

i = 1;
while i<length(tpCompMCpos)
    p1 = tpCompMCpos(i);
    p2 = tpCompMCpos(i+1);
    if (p2-p1)<=w
        if compMCtempMod(p1)>compMCtempMod(p2)
            tpCompMCpos(i+1)=0;
        else
            tpCompMCpos(i)=0;
        end
    end
    i=i+1;
end

tpCompMCpos = tpCompMCpos(tpCompMCpos>0);

clear i p1 p2

for i = 1:length(tpCompMCpos)-1
    p1 = tpCompMCpos(i);
    p2 = tpCompMCpos(i+1);
    
    if (p2-p1)>intvl && (p2-p1)<100
        if compMCtempMod(p1)>compMCtempMod(p2)
            tpCompMCpos(i+1)=0;
        end
    end
end

tpCompMCpos = tpCompMCpos(tpCompMCpos>0);

clear i

for i = 1:length(tpCompMCpos)-1
    
    p1 = tpCompMCpos(i);
    p2 = tpCompMCpos(i+1);
    
    if p2-p1<intvl && compMCtempMod(p1-1)==0
            tpCompMCpos(i)=0;
        
    end
end

tpCompMCpos = tpCompMCpos(tpCompMCpos>0);



for i = 1:length(tpCompMCpos)
    
    idk = [tpCompMCpos(i)-intvl:tpCompMCpos(i)+intvl];
    [c,ia,ib]=intersect(idk,tpCompMCpos);
    
    if length(c)~=1
        [maxval,maxid]= max(compMC(c));
        tpCompMCpos(i)=c(maxid);
    end
    
end

tp4v1 = [tpCompMCpos' 0];
tp4v2 = [0 tpCompMCpos'];

tp4v3 = tp4v1-tp4v2;
if ~isempty(tp4v3==0)
    tpCompMCpos(tp4v3==0)=0;
end


tpCompMCpos = tpCompMCpos(tpCompMCpos>0);
MC = tpCompMCpos;
cprintf('Hyperlinks','ValveTiming function: MC\n');
clear i
for i = 1:length(MC)
    temp = intensityCheck(MC(i),(sum(IP2)+sum(IN2)).*envP,'mc');
    
    if abs(diff([temp MC(i)]))<=2 || (temp-MC(i)<w && temp-MC(i)>0)
        if ~(s5v1(MC(i))>100 && s5v1(MC(i))>s5v1(temp))
            MC(i)=temp;
        else
            fprintf('orignial mc %d\n',MC(i));
        end
    elseif MC(i)-temp>2 && MC(i)-temp<w
        MC(i) = temp+round((MC(i)-temp)/2);
    elseif s3(temp)> s1(MC(i)) &&  (temp-MC(i)<=5 && temp-MC(i)>0)
        MC(i) = temp+round((MC(i)-temp)/2);
    end
    
    
end
clear i
MCorig = [MC';ones(1,length(MC))];

%%
if ~isempty(MC)
    wrongMCid =find(s5v2(MC)<max(s5v2(MC))/2);
    wrongMC = MC(wrongMCid);
    
    if ~isempty(wrongMC)
        
        for i = 1:length(wrongMC)
            mc = wrongMC(i);
            if ( mc-10>0 && mc+10<=N &&(mean(s1(mc-10:mc-5))>s1(mc) || s1(mc)<mean(s1(mc+5:mc+10))))...
                    || s5(mc)<1
                
                fprintf('%d might be wrong\n',wrongMC(i));
                
                blkid = mc-k1-w:mc+k1;
                blk = sumtempIP(blkid);
                [peaks,locs]=max(blk);
                
                if ~(abs(blkid(locs)- (mc+k1))<3 && (mc+k2<N && s1(mc+k1)/max(s1(mc+k1+w:mc+k2))<0.7))    
                    mc = blkid(locs)+1;
                end
                
                if s5v2(blkid(locs))>max(s5v2(MCorig(1,:)))/2
                    MC(wrongMCid(i))= mc;
                elseif mc-k1>0 && mc+k1<=N &&(mean(s1(mc-k1:mc-w))<s1(mc) && s1(mc)>mean(s1(mc+w:mc+k1)))
                    MC(wrongMCid(i))= mc;                           
                elseif mc-k1>0 && mc+w<=N &&(mean(s1(mc-k1:mc-w))<s1(mc) && s1(mc)>mean(s1(mc+1:mc+w)))
                    
                    MC(wrongMCid(i))= mc;
                else
                    MC(wrongMCid(i))= 0;
                    MCorig(2,wrongMCid(i))=0;
                end
            end
        end   
    end
    
end

MC=MC(MC>0);









