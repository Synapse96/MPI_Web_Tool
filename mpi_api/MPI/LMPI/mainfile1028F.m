%% This section needs to be run first if you want to run the following sections seperately
function [resultTable,resultValTime,MCnew,valueMC,AOnew,valueAO,ACnew,valueAC,MOnew,valueMO]=mainfile1028F(I_orig,time)
%%      
%{    
/** 
 * Define filters for :
 *     - image normalization 
 *     - envelope construction
**/    
%}

% -------------------------------------------------------------------------

d = designfilt('lowpassfir','DesignMethod','equiripple', ...
    'PassbandFrequency',0.05,'StopbandFrequency',0.1, ...
    'PassbandRipple',3,'StopbandAttenuation',60);

d2 = designfilt('lowpassfir','DesignMethod','equiripple', ...
    'PassbandFrequency',0.03,'StopbandFrequency',0.1, ...
    'PassbandRipple',3,'StopbandAttenuation',60);

d3 = designfilt('highpassfir','DesignMethod','equiripple', ...
    'PassbandFrequency',0.08,'StopbandFrequency',0.03, ...
    'PassbandRipple',3,'StopbandAttenuation',60);


%%      
%{    
/** 
 * Image filtering and envelope construction
**/    
%}

% -------------------------------------------------------------------------

    I = I_orig;
  
    I = imadjust(I);
    I = medfilt2(I,[5 5]);
    
    [counts,~]=imhist(I);
    [level1]=triangle_th(counts,256);
    
    [counts2,~]=imhist(I_orig);
    [level2]=triangle_th(counts2,256);
    
    clear  counts counts2 x x2 data1 data2
    
    [IP,IN,ipos,ineg,~,~,envP,envN] = constructEnvelope(I,level1);
    [IP1,IN1,~,~,~,~,envP1,~] = constructEnvelope(I_orig,level2);
    
    envN = movavgFilt(envN,5,'Center');
    envP = movavgFilt(envP,5,'Center');
    s5 = movavgFilt(sum(IP),5,'Center');
    s6 = movavgFilt(sum(IN),5,'Center');
    
    s1 = envP;
    s2 = envP.*ipos;
    s3 = envN;
    s4 = envN.*ineg;
    s6v2=sum(IN1);
    s1v2 = filtfilt(d,s1);
    
    N = size(I,2);
    T = 0.002;
    
    [Gmag1,Gdir1,Gx1,Gy1]= plotGmag(IP1,[],1:N,[],'no plot');
    [Gmag2,Gdir2,Gx2,Gy2]= plotGmag(IN1,[],1:N,[],'no plot');
 
    %%
    %{
    /**
     * Envelope signal segmentation
     * Valve event detection
    **/
    %}
    
    % -------------------------------------------------------------------------
    
    [AO,valAO,AC,valAC,MO,valMO,MC,valMC,s5new,sumtempIP,tempIN,MCorig,sumGx1,sumGx1Orig]= valveTimingsDetection(I_orig);
    
    clear  ipos maPos maNeg  ipos1 ineg1 maPos1 maNeg1  envN1
    
    f2 = s1-s3;
    f1 = s1+s3;
    f3 = f2; f3(f3<0)=0;
    f4 = f2; f4(f4>0)=0; f4 = abs(f4);
    
    f2v2 = filtfilt(d,f2);
    f1v2 = filtfilt(d,f1);
    
    f5 = s5+s6;
    wd = 0.4/T;
    
    [xRawPeakPoint,xRawBottomPoint,PeakDiffH,TroughDiffH, PulseWidth,PulseAmplitude]=PeakTrough(f2v2,wd);
    
    yRawPeakPoint = f2v2(xRawPeakPoint);
    yRawBottomPoint = f2v2(xRawBottomPoint);
    
    segmarker = xRawBottomPoint;
    if segmarker(1)> min(PulseWidth)/4*3 && segmarker(1)< max(PulseWidth)
        segmarker=[1 segmarker];
    elseif segmarker(1)> max(PulseWidth)
        segmarker=[segmarker(1)-max(PulseWidth)+1 segmarker];
        
    end
    
    if N-segmarker(end)> min(PulseWidth)/4*3 && N-segmarker(end)< max(PulseWidth)
        segmarker=[segmarker N];
    elseif N-segmarker(end)> max(PulseWidth)
        segmarker=[segmarker segmarker(end)+max(PulseWidth)];
    end
    
    clear k

%%      
%{    
/** 
 * Image segment quality classification
 *     - priamry: overall waveform quality
 *     - secondary: valve click quality
**/    
%}

% -------------------------------------------------------------------------

    k = 0;
    mpi = zeros(1,length(segmarker)-1);
    ict = zeros(1,length(segmarker)-1);
    irt = zeros(1,length(segmarker)-1);
    et = zeros(1,length(segmarker)-1);
    eventNo = zeros(1,length(segmarker)-1);
    mcintens = zeros(1,length(segmarker)-1);
    peakAODiff = zeros(1,length(segmarker)-1);
    aoVar= zeros(1,length(segmarker)-1);
    
    newfea1 = zeros(1,length(segmarker)-1);
    
    AOorig = AO;
    AOindict = AO(:,2);
    AO = AO(:,1);
    
    MOorig = MO;
    MOindict = MO(:,2);
    MO = MO(:,1);
    
    MC(s1(MC)==0)=[];
    AO(s1(AO)==0)=[];
    
    MCnew = zeros(1,length(segmarker)-1);
    AOnew = zeros(1,length(segmarker)-1);
    ACnew = zeros(1,length(segmarker)-1);
    MOnew = zeros(1,length(segmarker)-1);
    AOindictnew = nan(1,length(segmarker)-1);
    MOindictnew = nan(1,length(segmarker)-1);
    MCindictnew= nan(1,length(segmarker)-1);
    fprintf('\n');
    cprintf('*blue',  'Segment check\n');
    
    classCount = zeros(1,length(segmarker)-1);
    
    while k<length(segmarker)-1
        
        k=k+1;
        roi = segmarker(k):segmarker(k+1);
        mc = intersect(MC,roi); mcOrig = intersect(MCorig(1,:),roi);
        ao = intersect(AO,roi);
        ac = intersect(AC,roi);
        mo = intersect(MO,roi);moindict=[];
        topPeak = intersect(xRawPeakPoint,roi);
        fprintf('SEGMENT %d: \n',k);
        if length(ao)==1
            aoindict = AOorig(AOorig(:,1)==ao,2); AOindictnew(k) = aoindict;
        end
        
        %{
        /**
         * priamry:
         *     - threshold setup
         *     - overall waveform quality evaluation
        **/
        %}
        normalLength = 300;
        ejectMin = 65;
        valSepMin = 5;
        areaPerctgThreshold = 0.8;
        zeroAreaPerctgThreshold = 0.05;
        baselineVertical = 64;
        goodValveToWaveRatio = 0.6;
        waveLengthRatio = 0.5;
        definedBlock = 35;        

        lroi = length(roi);
        
        if lroi>normalLength
            fprintf('No good wave to get a proper cycle\n');
            classCount(k)=classCount(k)+inf;
            continue
        end
        
        roi1 = s1(roi);
        roi2 = s3(roi);
        zeroPerctg1=length(roi1(roi1==0))/lroi;
        zeroPerctg2=length(roi2(roi2==0))/lroi;
        cv1 = std(roi1)/mean(roi1);
        cv2 = std(roi2)/mean(roi2);
        
        smlPerctg1=length(roi1(roi1<10))/lroi;
        smlPerctg2=length(roi2(roi2<10))/lroi;
        
        lrgPerctg1=length(roi1(roi1>median(s1)/2))/lroi;
        lrgPerctg2=length(roi2(roi2>median(s3)/2))/lroi;
        
        if (smlPerctg1>=areaPerctgThreshold || smlPerctg2>=areaPerctgThreshold)
            fprintf('Upper or lower wave is missing\n');
            classCount(k)=classCount(k)+2;
            if length([mc' ao' ac' mo'])<=2
            MCnew(k)=0;
            AOnew(k)=0;
            ACnew(k)=0;
            MOnew(k)=0;
            mpi(k)=0;
            classCount(k)=classCount(k)+inf;
            continue
            end
        end
        
        if (lrgPerctg1>=areaPerctgThreshold || lrgPerctg2>=areaPerctgThreshold)
            fprintf('Nosiy baseline\n');
            classCount(k)=classCount(k)+2;
            if length([mc' ao' ac' mo'])<=2 && isempty(mc)&&isempty(ao)
            MCnew(k)=0;
            AOnew(k)=0;
            ACnew(k)=0;
            MOnew(k)=0;
            mpi(k)=0;
            classCount(k)=classCount(k)+inf;
            continue
            end
        end
        
        if ~(length(mc)==1&&length(ao)==1&&length(ac)==1&&length(mo)==1)
            classCount(k)=classCount(k)+1;
        end
        
        
        %{
        /**
         * priamry: overall waveform quality
         *     - clean up in edge roi
         *     - valve sequence order check in roi
         *     - solving misdetection caused by nosiy baseline
        **/
        %}
        
        if ~isempty(ac)&& abs(ac(1)-roi(1))<lroi/2
            ac(1) = [];
        end
        
        if length(mo)==1&& ...
                ((length([mc' ao' ac'])==3 && (length(mc)==1&&mo<mc))||(length([mc' ao' ac'])==2 && mo<min([mc' ao' ac'])))
            mo =[];
        end
        
        if  length(ao)==2&&length(mc)==1 &&length(mo)==1&&length(ac)==1&&(ao(1)<mc && ao(2)>mc)
            ao = ao(ac-ao>ejectMin & ac-ao<normalLength/3);
        end
        
        if length(ao)>1
            if ~isempty(ao((ao-roi(1))<lroi/4))
                ao((ao-roi(1))<lroi/4) = [];
            else
                xPeak = intersect(xRawPeakPoint,roi);
                ao(abs(xPeak - ao)>normalLength/10)=[];
            end
        end
        
        
        % solving misdetection caused by nosiy baseline     
        if length(mo)==2 || (length(mo)==1&&length(ac)==1&&abs(ac-mo)<valSepMin)
            if length(mc)==1 &&length(ao)==1&&length(ac)==1&&( ac-valSepMin*2>0 && ac+valSepMin*2<=N )
                if (max(s6(ac-valSepMin*2:ac-valSepMin))<s6(ac)&& s6(ac)>max(s6(ac+valSepMin:ac+valSepMin*2)))&&...
                        (max(s3(ac-valSepMin*2:ac-valSepMin))<s3(ac)&& s3(ac)>max(s3(ac+valSepMin:ac+valSepMin*2)))
                    
                    if length(mo)==1
                        mo=[];
                    elseif  length(mo)==2 && mc<ao
                        mo((abs(ac-mo)<valSepMin)|mo<mc)=[];
                    end
                end
            end
        end
        
        if length(mo)==1 && length(ac)==1 && abs(ac-mo)<valSepMin
            z=find(diff(f2v2(roi)>0)~=0)+1;
            z = roi(z);
            xPeak = intersect(xRawPeakPoint,roi);
            if ~isempty(xPeak)
                z = z(z>xPeak);
                if ac-z<=valSepMin*2
                    mo = [];
                end
            end
        end
        
        if length(mc)==1 && length(ao)==1&& abs(mc-ao)<valSepMin
            if  zeroPerctg1<zeroAreaPerctgThreshold && zeroPerctg2<zeroAreaPerctgThreshold
                fprintf('[%0.2f %0.2f] [%0.2f %0.2f] [%0.2f %0.2f] ', ...
                    zeroPerctg1,zeroPerctg2,cv1,cv2,mean(roi1),mean(roi2));
                
                mc = [];
                ao = [];
                MCnew(k)=0;
                AOnew(k)=0;
                fprintf('MC AO too close,too much baseline interfernce\n')
            end
        end
       
        
        %{    
        /** 
         * mitrial closure (mc) check
         *     - multiple mc detections
         *     - mc check based on ao when those two events are far away
        **/    
        %}   
        
        if length(mc)>1 && isempty(mc([s1v2(mc)./s1(mc)]>goodValveToWaveRatio)) && ~isempty(ao)       
           mc(mc>ao)=[];    
        end
        
        if length(mc)>1
            if ~isempty(mc(((roi(end)-mc)/lroi)<waveLengthRatio))
                if length(ao)==1 && abs(ao-(mc(((roi(end)-mc)/lroi)<waveLengthRatio)))<40
                    ao = [];
                end
                mc(((roi(end)-mc)/lroi)<waveLengthRatio)=[];
            elseif s3(roi(1))> baselineVertical/3
                mc(mc-roi(1)<50)=[];
            end
        end
        
        if length(mc)>1
            mc = [];
        end
        
        %{    
        /** 
         * solve mc and ao misdetection when
         *     - detected mc and ao are far away due to noise
        **/    
        %}          
        
        mcflag =0;aoflag=0;mcTemp = 0;
        if ~isempty(ao)&&~isempty(mc)
            if  (ao-mc)<- ejectMin/2 ||(ao-mc<=-valSepMin && s1(ao)<s1(mc))
                mc = []; mcflag = 1;
                
            elseif ao-mc>definedBlock && (ao+10<=N&& s1(ao)/max(s1(ao+valSepMin:ao+valSepMin*2))>0.3 && s1(ao)/max(s1(ao+valSepMin:ao+valSepMin*2))<0.65 )
                fprintf('MC %d is too far\n',mc);mc = [];
                
            elseif ao-mc>=definedBlock*(2/3) 
               
                if s1v2(mc)/s1(mc)<=0.5
                    if (s1(ao)/max(s1(ao+1:ao+15))>goodValveToWaveRatio&&s1(ao)/max(s1(ao+1:ao+15))<0.9) || s1(ao)/max(s1(ao+1:ao+15))<0.3
                        ao = [];
                    end
                    MCnew(k)= mc;
                else
                    mcTemp = mc;
                    
                    %{    
                    /** 
                     * evaluate mc quality
                    **/    
                    %}  
                    
                    if (mc-10>0&&mc+10<=N)&& (mean(s5(mc-10:mc-5))<s5(mc)&&mean(s5(mc+5:mc+10))<s5(mc))...
                            && (mean(s1(mc-10:mc-5))<s1(mc)&&mean(s1(mc+5:mc+10))<s1(mc))
                        if s3(mc)/max(s3(mc-10:mc-5))>0.65 && s3(mc)/max(s3(mc-10:mc-5))<1
                            mc = [];mcflag = 1;
                        elseif ao-mc>30 && mc-20>0 && (s3(mc)>max(s3(mc-20:mc-5))||(topPeak>ao&&topPeak-ao<20&&s1(topPeak)>baselineVertical/3)) 
                            mc = [];mcflag = 1;
                        else
                            MCnew(k)= mc;
                        end
                        
                    else
                        mc = [];mcflag = 1;
                    end
                end
                
            end
            

            
        end
       

        %% 
        %{    
        /** 
         * mc redetection 
         *   
        **/    
        %}  
        
        flag=0;
        if length(mc)==1 && length(ao)==1
            if abs(mc-ao) < valSepMin && min([mc ao])+definedBlock <= N && min([mc ao])-definedBlock>0
                cprintf('Errors' ,'AO MC close pixel %d; ',ao); ptemp = min([mc ao]);
                blkid = min([mc ao])-definedBlock:min([mc ao])+definedBlock;
                blk = sumGx1(blkid);
                
                [peaks,locs]=findpeaks(blk);
                locs = blkid(locs);
                locs(peaks<2*median(blk(blk>0)))=[];
                
                [minval,minid]=min(abs(locs-ptemp));
                if length(locs)>1 && minid == length(locs)
                    mc = [];flag = 1; cprintf('Errors' ,' this is AO ');
                    
                elseif length(locs)>1 && minid == 1
                    ao = [];cprintf('Errors' ,' this is MC ');
                else
                    cprintf('Errors' ,'Cannot distinguish ');
                end
                fprintf('\n');
            end
        end
        
        if length(ao)==1 && isempty(mc)
            fprintf('MC missing\n');
            AOnew(k) = ao; 
            classCount(k)=classCount(k)+1;

            s1flag=0;
            if ao-(definedBlock-valSepMin)>0
                aoid = ao-(definedBlock-valSepMin);
            else aoid = 1;
            end
            
            aoid = aoid:ao-valSepMin;
            widP = sumGx1Orig(aoid);
            [mcval,mcP] = findpeaks(widP);
            mcP(mcval<max(mcval)/3)=[];mcval(mcval<max(mcval)/3)=[];
            [ms1,ms1id]= findpeaks(s1(aoid));
            
            %{    
             /** 
              * MCindictnew is the quality indicator for classification purpose
              *   
             **/    
            %} 
            
            if length(mcP)==1&&length(ms1id)==1&&abs(ms1id - mcP)<=3
                fprintf('MC is empty: only 1 good local peak\n');
                mc = mcP+1;
            else
                if length(mcP)>=3
                    fprintf('MC is empty: use ao to fine mc, but the intensity is not very obvious, %d local peaks\n',length(mcP));
                    MCindictnew(k) = 3;
                end
                MCindictnew(k) = 2;
                mcPorig = mcP;
                mcP = mcP(mcval==max(mcval));
                mcval = max(mcval);
             
                if isempty(ms1id) 
                    mc = mcP+1;MCindictnew(k) = 2;
                elseif sum(abs(ms1id - mcP)<=3)~=0
                    if length(ms1id)>1 || mcP>ms1id
                        mc = mcP+1;MCindictnew(k) = 2;
                    else
                        mc = ms1id(abs(ms1id - mcP)<=3);s1flag = 1;MCindictnew(k) = 2;
                    end
                else
                    fprintf('MC is empty: cannot use upper wave, consider lower wave\n');MCindictnew(k) = 3;
                    widN = tempIN(aoid);
                    [mcvalN,mcN] = findpeaks(widN);
                    mcN(mcvalN<max(mcvalN)/3|mcvalN<1)=[];
                    mcvalN(mcvalN<max(mcvalN)/3|mcvalN<1)=[];          
                    if  sum(abs(mcN - mcP)<=2)~=0
                        mc = mcP+1;
                    elseif length(ms1id)==1&& sum(abs(ms1id - mcN)<=3)~=0
                        mc = ms1id;s1flag = 1;
                    end
                end     
            end

            if ~isempty(mc)
                mc = aoid(mc);
                if abs(intensityCheck(mc,s5)-mc)<=2&&s1flag==0
                    mc = intensityCheck(mc,s5);
                end
            else
                cprintf('Error','cannot fine mc\n');
            end
            
            
            if ~isempty(mc)
                if mc<=mcTemp && mcflag ==1
                    mc = [];MCnew(k)=0;
                elseif length(mc)==1&& ao-mc>valSepMin&& ao-mc<definedBlock - valSepMin 
                    MCnew(k)= mc;
                else
                    if length(mc)==1&& (ao-mc<=valSepMin|| ao-mc>=definedBlock - valSepMin )
                        cprintf('Error','redetected mc too close or too far\n');
                    end
                    mc = [];MCnew(k) = 0;
                end
            end
        end
        
        
        %{
         /**
          * Use original indicator MCorig(2,:)
          *
         **/
        %}
        
        if length(mcOrig)==1 && length(mc)==1&& flag~=1
            if abs(diff([mc,mcOrig]))<=2 && MCorig(2,MCorig(1,:)==mcOrig)==0
                aoid = mc:ao-5;
                wid = sumtempIP(aoid);
                [mcval,mc] = findpeaks(wid);
                mc(mcval<max(mcval)/3)=[];mcval(mcval<max(mcval)/3)=[];
                mc(mcval<0.5)=[];
                if ~isempty(mc)
                    mc = aoid(mc(mcval==max(mcval)))+1;
                else mc = mcOrig;
                end
                
            elseif s3(mcOrig)<s3(mc)&& s1(mcOrig)>s1(mc)&& mcflag==0 
                mc = mcOrig;
            end
        end
        
        if ~isempty(mc)
            if (~isempty(ao)&&ao-mc<10) && (mc-5>0&&mean(s3(mc-5:mc))<1)
                fprintf('MC: detected mc %d too close to ao, not right, cannot find a proper mc\n',mc);
                mc = 0;MCnew(k)=[];
            elseif length(mcOrig)==1 && length(mc)==1&&abs(diff([mc,mcOrig]))<=2 && mcflag ==1
                fprintf('MC: detected mc %d too far from ao, not right\n',mc);
                mc = 0;MCnew(k)=[];
            end
        end
        
        
        acOrig=[];
        if length(ac)==1 && ~isempty(ao)
            if ac-ao < ejectMin || ac-ao>normalLength/3
                acOrig=ac;
                ac = [];
            end
        end
        
        %{
         /**
          * segment contains another mitral opening
          *
         **/
        %}
        
        if k==1 && length(mo)>1  
            if ~isempty([mc' ao' ac'])
                mo = mo(mo>min([mc' ao' ac']));
            end
        elseif length(mo)>1 && (length(mc)==1&&length(ao)==1&&length(ac)==1)
            pause
            if (ao-mc>valSepMin&&ao-mc<definedBlock)&&(ac-ao>ejectMin&&ac-ao<normalLength/3)
                mo = mo(mo-ac<35);
            end
        end
        
        
        if isempty(ac) && isempty(mo)
           classCount(k)=classCount(k)+1;

            fprintf('No AC or MO');
            fprintf('(ao is %d); ',ao);
            if ~isempty(ao)
                region = ao+round(length(ao:roi(end))/3): roi(end);
            else
                region = roi(round(lroi/3):end);
            end
            ac = findAC(s4,d3,region);
            
            if abs(ac - region(end))<10 || (length(ao)==1&& ac-ao>normalLength/3)
                ac = findAC(s3,d3,region(1):region(end)-10);
                if (tempIN(ac)<0.5&&tempIN(ac-1)<0.5) || (length(ao)==1&&ac-ao<ejectMin)
                    ac = findAC(s3+s1,d3,region);
                end
            elseif abs(ac - region(1))<10
                ac = findAC(s3,d3,region(10):region(end));
            end
            
            
            if (length(ao)==1&& (ac-ao>normalLength/3 || ac-ao<ejectMin)) || (length(mc)==1&& ac-mc>normalLength/3 + valSepMin*4 )
                ac = [];
            else
                if (isempty(ao)|| isempty(mc)) && (~(isempty(ao)&&isempty(mc)))% if ao or mc is missing, won't go into the loop
                    if ac+definedBlock>N
                        widEnd = N;
                    else widEnd = ac+definedBlock;
                    end
                    
                    wid = s6(ac+valSepMin:widEnd);
                    [moval,mo] = max(wid);
                    
                    mo = mo+ac+valSepMin-1;
                    blkid = (ac+valSepMin*2):mo;
                    blknew = s6(blkid);
                    difblk = gradient(blknew);
                    [maxval,maxid]=max(difblk);
                    mo=blkid(maxid)+1;
                    [mo,moindict]=refineAndClassifyMO(tempIN,mo,s3,N,s6v2);mo = mo+1;
                end
            end
            
            if lrgPerctg2>=areaPerctgThreshold && ~isempty(ac)
                if ~(ac-valSepMin>0&& ac+valSepMin<=N && (mean(s3(ac-valSepMin:ac-1))<s3(ac)&&  mean(s3(ac+1:ac+valSepMin))<s3(ac))&& ...
                        (mean(s6(ac-valSepMin:ac-1))<s6(ac)&&  mean(s6(ac+1:ac+valSepMin))<s6(ac)))
                    fprintf('AC: too much large magnitude baseline noise, not sure if AC %d is right',ac);
                classCount(k)=classCount(k)+inf;
    
                end
            end
            fprintf('\n');
        end
        
        
        if length(mo)==1 && length(ac)==1 && mo-ac>definedBlock
            if mo+10<=N && s3(mo)/max(s3(mo+valSepMin:mo+valSepMin*2))>0.4
                
                widID = mo-definedBlock:mo-valSepMin;
                wid = tempIN(mo-definedBlock:mo-valSepMin);
                [acval,ac] = findpeaks(wid);
                ac(acval<max(acval)/3)=[];
                acval(acval<max(acval)/3)=[];
                ac = ac+mo-definedBlock-1;
                
                widP = sumGx1Orig(widID);
                [acvalP,acP] = findpeaks(widP);
                acP(acvalP<max(acvalP)/3)=[];
                acvalP(acvalP<max(acvalP)/3)=[];
                acP(acvalP<1)=[];
                acvalP(acvalP<1)=[];
                if ~isempty(acP)
                    acP = widID(acP);
                end
                
                [acvalPN, acPN] = max(wid+widP); acPN = widID(acPN);
                if sum(abs(acPN-ac)<3)~=0 || (~isempty(acP)&& sum(abs(acPN-acP)<3)~=0 )
                    ac = acPN+1;
                end
            end
        elseif length(mo)==1 && length(ac)==1 && mo<ac
            if s3(roi(end))>20 && roi(end)-mo>50
                mo = [];moindict = [];
            end
        end
        if length(mo)==1
            if isempty(moindict)
            moindict = MOorig(MOorig(:,1)==mo,2);
            end
            MOindictnew(k)=moindict;
        end

        %%
        %{    
        /** 
         * after each of the valve click is checked 
         *     - generate class quality score
         *     - put the checked valve event into corresponding sequence 
        **/    
        %}  
        
        if length([mc' ao' ac' mo'])<=1
            
            MCnew(k) =  0;
            AOnew(k) =  0;
            ACnew(k) =  0;
            MOnew(k) =  0;
            classCount(k)=classCount(k)+inf;

        elseif length(mc)==1 && length(ao)==1 && length(ac)==1 && length(mo)==1
            ict(k) = ao-mc;
            irt(k) = mo-ac;
            et(k) = ac-ao;
            
            mpi(k) = (ao-mc+mo-ac)/(ac-ao);
            MCnew(k) =  mc;
            AOnew(k) =  ao;
            ACnew(k) =  ac;
            MOnew(k) =  mo;
        elseif length(ao)==1 && length(ac)==1 && length(mo)==1 && isempty(mc)
            
            MCnew(k) =  0;
            AOnew(k) =  ao;
            ACnew(k) =  ac;
            MOnew(k) =  mo;
            classCount(k)=classCount(k)+inf;
            
        elseif length(mc)==1 && length(ac)==1 && length(mo)==1 && isempty(ao)
            classCount(k)=classCount(k)+1;
            
            fprintf('AO missing for %d\n', mc)
                widID = mc+valSepMin:mc+25;         
                [ao,aoindict] = refineAndClassifyAO(widID,sumGx1,sumGx1Orig,N,mc+25,s1,envP1);
                
                MCnew(k) =  mc;
                ACnew(k) =  ac;
                MOnew(k) =  mo;
                
                if ao~=0
                    AOnew(k)= ao;
                    ict(k) = ao-mc;
                    irt(k) = mo-ac;
                    et(k) = ac-ao;
                    mpi(k) = (ao-mc+mo-ac)/(ac-ao);
                    
                    AOindictnew(k) = aoindict;
                else
                    ao = [];AOnew(k)= [];
                    classCount(k)=classCount(k)+inf;
                    AOindictnew(k) =[];
                end
            
        elseif length(ao)==1 && length(ac)==1 && length(mc)==1 && isempty(mo)
            classCount(k)=classCount(k)+1;
            if ac+definedBlock>N
                widEnd = N;
            else widEnd = ac+definedBlock;
            end
            
            wid = s6(ac+valSepMin:widEnd);
            [moval,mo] = max(wid);
            
            mo = mo+ac+valSepMin-1;
            blkid = (ac+valSepMin*2):mo;
            blknew = s6(blkid);
            difblk = gradient(blknew);
            [maxval,maxid]=max(difblk);
            mo=blkid(maxid)+1;moOrig = mo;
            
            [mo,moindict]=refineAndClassifyMO(tempIN,mo,s3,N,s6v2);
            mo=mo+1;
            
            if isempty(mo) || length(mo)>1
                mo = 0;moindict = nan;
                irt(k)=0;
                mpi(k)=0;
                classCount(k)=classCount(k)+inf;
            else
                
                ict(k) = ao-mc;
                irt(k) = mo-ac;
                et(k) = ac-ao;
                mpi(k) = (ao-mc+mo-ac)/(ac-ao);
            end
            
            MCnew(k) =  mc;
            AOnew(k) =  ao;
            ACnew(k) =  ac;
            MOnew(k) =  mo;
            MOindictnew(k)=moindict;
        elseif (length(mc)==1 && length(ao)==1 && length(mo)==1 && isempty(ac)) ||(length(mo)==1&&length(MO)>=length(segmarker)&&isempty(ac))
            classCount(k)=classCount(k)+1;
            if mo-definedBlock>0
                widID = mo-definedBlock:mo-valSepMin;
                wid = tempIN(mo-definedBlock:mo-valSepMin);
                [acval,ac] = findpeaks(wid);
                ac(acval<max(acval)/3)=[];
                acval(acval<max(acval)/3)=[];
                ac = ac+mo-definedBlock-1;
                
                wid2 = s6v2(widID);
                [acval2,ac2] = max(wid2);
                ac2 = widID(ac2);
                
                %{    
                 /** 
                  * For weak ac detection and validation  
                  *     - note that the original closure click can be weak
                  *     - closure clicks normally have reflections
                  *     - use both original and relection for better re
                  detection
                  **/    
                %}  
                widP = sumGx1Orig(widID);
                [acvalP,acP] = findpeaks(widP); acP = widID(acP);
                acP(acvalP<max(acvalP)/3)=[];
                acvalP(acvalP<max(acvalP)/3)=[];
                
                
                acP(acvalP<1)=[];
                acvalP(acvalP<1)=[];
            
                [acvalPN, acPN] = max(wid+widP); acPN = widID(acPN);
                
                if s6v2(acPN)<0.5 && mean(s1(acPN-10:acPN-5))>s1(acPN)
                    ac = ac2;
                else
                    
                    if sum(abs(acPN-ac)<3)~=0 || (~isempty(acP)&& sum(abs(acPN-acP)<3)~=0)
                        if ~isempty(acOrig)&&abs(acPN+1-acOrig)<3
                            ac = acOrig;
                        else
                            if s1(ac2)>baselineVertical/3 && ((ac2-acPN)>0&&(ac2-acPN)<4)
                                ac = acPN+2;
                            else
                                ac = acPN+1;
                            end
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isempty(ao)
                    ao = nan;classCount(k)=classCount(k)+inf;
                end
                if isempty(mc)
                    mc = nan;classCount(k)=classCount(k)+inf;
                end
                if length(ac)>1
                    if length(acP)==1&& sum(acP==ac)==1
                        ac = acP;
                    else
                        ac = nan;classCount(k)=classCount(k)+inf;
                    end
                end
                
                ict(k) = ao-mc;
                irt(k) = mo-ac;
                et(k) = ac-ao;
                mpi(k) = (ao-mc+mo-ac)/(ac-ao);
                
                MCnew(k) =  mc;
                AOnew(k) =  ao;
                ACnew(k) =  ac;
                MOnew(k) =  mo;
            end
            
            
        elseif length(mc)==1 && length(ao)==1 && isempty(ac) && isempty(mo)
            MCnew(k) =  mc;
            AOnew(k) =  ao;classCount(k)=classCount(k)+inf;
        else
            if isempty(ac)|| length(ac)~=1
                ac = 0;classCount(k)=classCount(k)+inf;
            end
            
            if isempty(mo) || length(mo)~=1
                mo = 0;classCount(k)=classCount(k)+inf;
            end
            ACnew(k) =  ac;
            MOnew(k) =  mo;
            
        end
        
        if MCnew(k)~=0 && (~isnan(MCnew(k)))
            mcintens(k) = s5(MCnew(k));
        else
            mcintens(k)=0;
        end
        eventNo(k) = length(mc)+length(ao)+length(ac)+length(mo);
        
        if ~isempty(ao)&& length(ao)==1 && (~isnan(ao))
            newfea1(k)= s1(ao)/(max(s1(ao:roi(end))));            
            if ~isempty(intersect(xRawPeakPoint,roi))
                aoVar(k)= s1(ao)/s1(intersect(xRawPeakPoint,roi));
            end
        end
        
        fprintf('\n');
   
    end
    
    
    %%
    %{
    /**
     * Generate new images for display purpose
     *     - transform the raw image to ultrasound like images
     *     - add detected valve events as vertical lines
    **/
    %}
    
    % -------------------------------------------------------------------------
    
    MCnew(isnan(MCnew))=0;
    MOnew(isnan(MOnew))=0;
    AOnew(isnan(AOnew))=0;
    ACnew(isnan(ACnew))=0;
    
    MCnew = MCnew(MCnew~=0);
    MOnew = MOnew(MOnew~=0);
    AOnew = AOnew(AOnew~=0);
    ACnew = ACnew(ACnew~=0);
    
    
    ultra = replicateUltrasound(I_orig);
    ultra = imadjust(ultra,[0.05 0.7],[0 0.69],1.2);
    
    mcOrig = MCnew;
    cprintf('*Blue','PostCheck MC\n');
    [mcNew,MCindict] = postCheckMC(MCnew,IP1,IN1,Gx1,Gx2,s5new,s1,AOnew);
    
    mcNewId = find(mcOrig~=mcNew);
    
    if ~isempty(mcNewId)
        k =1;
        while k<length(segmarker)
            
            roi = segmarker(k):segmarker(k+1);
            mc = intersect(mcNew,roi);
            ao = intersect(AOnew,roi); 
            ac = intersect(ACnew,roi);
            mo = intersect(MOnew,roi);
            mpinew = (ao-mc+mo-ac)/(ac-ao);
            numstr = sprintf('%0.2f ',mpinew);
            cprintf('Magenta', numstr);
            if isempty(mpinew)
                mpinew = 0;
            end
            mpi(k)=mpinew;
            k = k+1;
        end
        fprintf('\n')
    end
    MCindict = MCindict(mcNew>0);
    MCNew = mcNew(mcNew>0);
    
    plotLine = 1:2:floor(size(ultra,1)/2)*2;
    if ~isempty([AOnew MOnew])
        mask = ultra<0;
        mask(plotLine,[AOnew MOnew]) = 1;
        addOP = imoverlay(ultra,mask,[0 1 0]);
    else addOP = ultra;
        
    end
    if ~isempty([MCnew ACnew])
        mask = ultra<0;
        mask(plotLine,[MCnew ACnew]) = 1;
        addOPED = imoverlay(addOP,mask,[1 1 0]);
    else   addOPED = addOP ;
    end
    
    
    labels = {'ClassCount',[classCount]; 'AOindict',[AOindictnew];'MCindict',[MCindict];'MCindictnew',[MCindictnew];'MOindict',[MOindictnew]};
    
    valueAO = envP(AOnew);
    valueAC = envN(ACnew);
    valueMO = envN(MOnew);
    valueMC = envP(MCnew);
    
    %%
    %{
    /**
     * Put final results into the table, including:
     *     - valve event (pixel index)
     *     - time interval (ms)
     *     - MPI value
    **/
    %}
    
    % -------------------------------------------------------------------------
    
    clear k
    resultNum = nan(length(segmarker),4);
    resultValTime = nan(length(segmarker),4);
    for k = 1:length(segmarker)-1
        seg = segmarker(k):segmarker(k+1);
        
        mc = intersect(mcNew,seg);
        if isempty(mc)
            mc =0; mcR = 0;
        else mcR = time(mc);
        end
        
        ao = intersect(AOnew,seg);
        if isempty(ao)
            ao =0; aoR = 0;
        else aoR = time(ao);
        end
        
        ac = intersect(ACnew,seg);
        if isempty(ac)
            ac =0; acR = 0;
        else acR = time(ac);
        end
        mo = intersect(MOnew,seg);
        if isempty(mo)
            mo =0; moR = 0;
        else moR = time(mo);
        end
        
        resultNum(k,:) = [mc ao ac mo];
        resultValTime(k,:) = [mcR aoR acR moR];
    end
    
    resultNum(resultNum==0)=NaN;
    resultValTime(resultValTime==0)=NaN;
    
    resultNum(:,5) = (resultValTime(:,2)-resultValTime(:,1))*1000;
    resultNum(:,6) = (resultValTime(:,3)-resultValTime(:,2))*1000;
    resultNum(:,7) = (resultValTime(:,4)-resultValTime(:,3))*1000;
    resultNum(:,8) = (resultNum(:,5)+resultNum(:,7))./resultNum(:,6);
    
    hr = diff(resultNum(:,1));
    if length(MC)>length(MCnew)&& MC(end)-MCnew(end)>normalLength/2 && MC(end)-MCnew(end)<normalLength
        if isnan(hr(end))
            resultNum(1:length(hr),9) = [hr(1:end-1); MC(end)-MCnew(end)];
        end
    else
        resultNum(1:length(hr),9) = hr;
    end
    
    resultNum = num2cell(resultNum);
    labelTit = {'MC' 'AO' 'AC' 'MO' 'ICT' 'ET' 'IRT' 'MPI' 'HR'};
    resultTable = [labelTit;resultNum];
end
