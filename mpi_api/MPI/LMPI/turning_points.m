function [confirmed]=turning_points(x,threshold)
%Original version Stephen Redmond
%
%Modified by Philip de Chazal 30/5/07

x=x(:); %must be a column vector

x=[x(:);x(end)+eps;x(end)];

% finds the locations of peaks and troughs
% returns 

%tps=1 when at a peak and tps=-1 when at a through, tps=0 elsewhere
tps=[0;-sign(diff(sign(diff(x))));0];
tpidx=find(tps~=0);

%index of all turning points
pkth=tps(tpidx); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start searching for turning point using threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
inpeak=0;
possibleidx=[];
ref=x(1);
confirmed=NaN*ones(length(pkth),1);
k=1;
while i<=length(tpidx)

    %The aim of the following code is to eliminate all the local peaks and
    %troughs. A local peak or trough occurs when the height difference
    %between a peak and trough is less than 'threshold's
    
    %find first pk/tr
    if(inpeak==0 & abs(x(tpidx(i))-ref)>threshold)
        inpeak=pkth(i);
        possibleidx=tpidx(i);
		ref=x(tpidx(i));
    end
    
    %if looking for peak
    if( inpeak==1 & (ref-x(tpidx(i)))>threshold )
        %peak found when next trough is more then threshold away from
        %current peak
        confirmed(k)=possibleidx*tps(possibleidx);
        k=k+1;
        %this lower point could be next trough
        possibleidx=tpidx(i);
        ref=x(tpidx(i));
        inpeak=-1;
    elseif( inpeak==1 & x(tpidx(i)) > ref ) 
        possibleidx=tpidx(i);
        ref=x(tpidx(i));
    end
    
    %if looking for trough
    if( inpeak==-1 & (x(tpidx(i))-ref)>threshold )
        %trough found when next peak is more then threshold away from
        %current trough
        confirmed(k)=possibleidx*tps(possibleidx);
        k=k+1;
        %this higher point could be next peak
        possibleidx=tpidx(i);
        ref=x(tpidx(i));
        inpeak=1;
    elseif( inpeak==-1 & x(tpidx(i)) < x(possibleidx) ) 
        possibleidx=tpidx(i);
        ref=x(tpidx(i));
    end
       
    

    
%     if(~isempty(confirmed))
%         plot(idx,x(idx));
%         hold on;
%         idx2=find(abs(confirmed)>=idx(1) & abs(confirmed)<=idx(end));
%         plot(abs(confirmed(idx2)),x(abs(confirmed(idx2))),'g*')
%         plot(possibleidx,x(possibleidx),'black*')
%         plot(tpidx(i),x(tpidx(i)),'r*')
%         hold off;
%         pause;
% 	end
    i=i+1;
end
confirmed=confirmed( (~isnan(confirmed)));