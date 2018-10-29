function [xRawPeakPoint,xRawBottomPoint,PeakDiffH,TroughDiffH, PulseWidth,PulseAmplitude]=PeakTrough(signal,wd)

top=sortfilt1(signal,wd,100);
bottom=sortfilt1(signal,wd,0);

% to find the peaks
k=1;
for m=1:length(signal)
    if top(m)==signal(m)
        xRawPeakPoint(k)=m;
        yRawPeakPoint(k)=top(m);
        k=k+1;
    end
end

% to allocate 2 troughs for every peak
s=1;
for m=1:(length(xRawPeakPoint)-1)
    MinVal=min(signal(xRawPeakPoint(m):xRawPeakPoint(m+1)));
    for c=xRawPeakPoint(m):xRawPeakPoint(m+1)
        if signal(c)==MinVal
            xRawBottomPoint(s)=c;
            yRawBottomPoint(s)=MinVal;
            s=s+1;
        end
    end
end

xPeakPoint=xRawPeakPoint(2:length(xRawPeakPoint)-1);
yPeakPoint=yRawPeakPoint(2:length(yRawPeakPoint)-1);


for m=1:length(xPeakPoint)
    xBottomPoint(2*m-1)=xRawBottomPoint(m);
    yBottomPoint(2*m-1)=yRawBottomPoint(m);
    xBottomPoint(2*m)=xRawBottomPoint(m+1);
    yBottomPoint(2*m)=yRawBottomPoint(m+1);
end

xPeakPoint11=xPeakPoint;
xBottomPoint11=xBottomPoint;

TroughDiffH = zeros(1,length(xPeakPoint));
PulseWidth = zeros(1,length(xPeakPoint));
PulseAmplitude = zeros(1,length(xPeakPoint));

for m=1:length(xPeakPoint)
    TroughDiffH(m)=abs(yBottomPoint(2*m)-yBottomPoint(2*m-1));
    PulseWidth(m)=abs(xBottomPoint(2*m)-xBottomPoint(2*m-1));
    PulseAmplitude(m)=max([abs(yPeakPoint(m)-yBottomPoint(2*m-1)) abs(yPeakPoint(m)-yBottomPoint(2*m))]);
end

clear m

PeakDiffH= zeros(1,length(xRawPeakPoint)-1);

for m = 1:length(xRawPeakPoint)-1
    PeakDiffH(m) = yRawPeakPoint(m+1)-yRawPeakPoint(m);  
end



