function [MC,MCindict] = postCheckMC (MCnew,IP1,IN1,Gx1,Gx2,s5new,s1,AOnew)
% Post checking the quality of detected MC clicks, based on the
% shape/intensity 

MC = MCnew;
flagClose = zeros(1,length(MC));
MCindict = NaN(1,length(MC));
w = 5;
k = 10;

if ~isempty(MCnew)
    
    clear m pmc
    for m = 1:length(MCnew)
        
        pmc = MCnew(m);
        blkPmcID = pmc-5:pmc+5;
        blkPmcIP1 = IP1(:,blkPmcID);
        blkPmcIN1 = IN1(:,blkPmcID);
        blkPmcGx1 = Gx1(:,blkPmcID);
        blkPmcGx2 = Gx2(:,blkPmcID);
        [maxblkPmcIP1,maxlocblkIP1] = max( sum(blkPmcIP1));
        [maxblkPmcGx1,maxlocblkGx1] = max( sum(blkPmcGx1));   
        
        blkPmcGx1( blkPmcGx1<0)=0;
        blkPmcGx2(blkPmcGx2<0)=0;
        
        pmcInten = s5new(pmc);
        pmcBefore = mean(s5new(pmc-5:pmc-1));
        pmcAfter = mean(s5new(pmc+1:pmc+5));
        
        if pmc-k>0&&pmc+k<=length(s1)
            pmcBefore2 = mean(s5new(pmc-k:pmc-w));
            pmcAfter2 = mean(s5new(pmc+w:pmc+k));
        end
        
        pmcEnvlp = s1(pmc);
        pmcEnvlpBefore = mean(s1(pmc-w:pmc-1));
        pmcEnvlpAfter = mean(s1(pmc+1:pmc+w));
        
        fprintf('[%0.2f %0.2f %0.2f]\n[%0.2f %0.2f %0.2f] [%d %d] [%0.2f ] [%0.2f]',...
            pmcEnvlpBefore,pmcEnvlp, pmcEnvlpAfter,...
            pmcBefore, pmcInten,pmcAfter,...
            maxlocblkIP1,maxlocblkGx1,...
            std(sum(blkPmcIP1))/mean(sum(blkPmcIP1)), ...
            std(s1(pmc-w:pmc+w))/mean(s1(pmc-w:pmc+w)));
       
        sig = sum(blkPmcGx1);
        id = pmc-w:pmc+w;
        [peaks,locs] = findpeaks(sig);
        
        if (pmcEnvlpBefore<=pmcEnvlp && pmcEnvlp>=pmcEnvlpAfter) &&...
                (pmcBefore<pmcInten && pmcInten>pmcAfter)
            if pmc-k>0&&pmc+k<=length(s1)
                cprintf('Error','\n[%0.2f %0.2f %0.2f]',pmcBefore2, pmcInten,pmcAfter2);
                MCindict(m)=0;
            end
            fprintf('\n%d seems to be a peak in both envelop and intensity\n',pmc);
            
        elseif pmcEnvlpBefore>pmcEnvlp && pmcEnvlp<pmcEnvlpAfter
            fprintf('\n%d envelope smaller than both before and after areas\n',pmc)
            if pmcBefore<pmcInten && pmcInten>pmcAfter
                if pmc-k>0&&pmc+k<=length(s1)
                    cprintf('Error','[%0.2f %0.2f %0.2f] ',pmcBefore2, pmcInten,pmcAfter2);
                    if  pmcBefore2<pmcInten/3 &&  pmcAfter2<pmcInten/3*2
                        MCindict(m)=2;
                    else MCindict(m)=3;
                    end
                end
                fprintf('but the intensity is bigger than both areas\n');
            else
                if length(locs(peaks>max(peaks)/2))==1
                    pmctemp = id(locs(peaks>max(peaks)/2)); MCindict(m)=3;
                    if s5new(pmctemp)> mean(s5new(pmctemp-w:pmctemp-1))&&...
                            s5new(pmctemp)> mean(s5new(pmctemp+1:pmctemp+w))
                        MC(m)=pmc+round(diff([pmc pmctemp])/2);
                        fprintf('changed MC to %d according to gradient 62\n',MC(m));
                    end
                else fprintf('even intensity is small, might be wrong\n');
                    MCindict(m)=inf;
                end
            end
            
        elseif pmcEnvlpBefore>pmcEnvlp && pmcEnvlp>pmcEnvlpAfter
            fprintf('\n%d seems to be together with bright junk, envelope before is bigger\n',pmc);
            MCindict(m)=2;
            if pmcBefore<pmcInten && pmcInten>pmcAfter
                if pmc-10>0&&pmc+10<=length(s1)
                    cprintf('Error','[%0.2f %0.2f %0.2f] ',pmcBefore2, pmcInten,pmcAfter2);
                    if  pmcBefore2<pmcInten/3 &&  pmcAfter2<pmcInten/3*2
                        MCindict(m)=2;
                    else MCindict(m)=3;
                    end
                end
                fprintf('but the intensity is bigger than both areas\n');
            else
                if length(locs(peaks>max(peaks)/2))==1
                    pmctemp = id(locs(peaks>max(peaks)/2));
                    if abs(diff([pmctemp pmc]))<=2
                        fprintf('might be right, not far from gradient peak\n');MCindict(m)=3;
                    end
                else
                    fprintf('even intensity is small, might be wrong\n'); MCindict(m)=inf;
                end
            end
            
        elseif  pmcEnvlpBefore<pmcEnvlp && pmcEnvlp<pmcEnvlpAfter
            fprintf('\n%d envelop bigger than before, smaller than after\n',pmc); flagClose(m) = 1;
            
            if sum(abs(AOnew-pmc)<10)==1
                fprintf('too close to following AO\n');
                MC(m)=0; MCindict(m)=inf;
            else
                
                if pmcBefore<pmcInten && pmcInten>pmcAfter
                    if pmc-10>0&&pmc+10<=length(s1)
                        cprintf('Error','[%0.2f %0.2f %0.2f] ',pmcBefore2, pmcInten,pmcAfter2);
                        if  pmcBefore2<pmcInten/3 &&  pmcAfter2<pmcInten/3*2
                            MCindict(m)=2;
                        else MCindict(m)=3;
                        end
                    end
                    fprintf('but the intensity is bigger than both areas\n');
                    
                else
                    
                    if length(locs(peaks>max(peaks)/2))==1
                        pmctemp = id(locs(peaks>max(peaks)/2));MCindict(m)=3;
                        if s5new(pmctemp)> mean(s5new(pmctemp-w:pmctemp-1))&&...
                                s5new(pmctemp)> mean(s5new(pmctemp+1:pmctemp+w))
                            MC(m)=pmctemp;
                            fprintf('changed MC to %d according to gradient \n',pmctemp);
                        end
                    else  fprintf('might be wrong\n');MC(m)=0;MCindict(m)=inf;
                    end
                end
            end
            
            
        elseif  (pmcEnvlpBefore>pmcEnvlp || pmcEnvlp<pmcEnvlpAfter) && ...
                ( pmcBefore > pmcInten ||  pmcInten<pmcAfter)
            
            if round(abs(pmcInten-pmcAfter))<=1 ||  round(abs(pmcEnvlp-pmcEnvlpAfter))<=1
                if  abs(diff([maxlocblkIP1 maxlocblkGx1]))>w
                    fprintf('\n%d not correct\n?',pmc);MCindict(m)=inf;
                else
                    fprintf('\n%d might be correct?\n',pmc);MCindict(m)=4;
                end
            else
                fprintf ('\nnot correct?\n');MCindict(m)=inf;
            end
            
        elseif (pmcBefore<pmcInten && pmcInten>pmcAfter) && ...
                (round(pmcInten-pmcBefore)<=1 || round(pmcInten-pmcAfter)<=1)
            
            if round(pmcInten-pmcBefore)<=0.5 && round(pmcInten-pmcAfter)<=0.5
                fprintf('\nIntensity change not obvious, not correct?\n')
                MCindict(m)=inf;
            else
                fprintf('\nnot sure?');
                fprintf('\n%d %d',std(sum(blkPmcIP1))/mean(sum(blkPmcIP1)),(pmcInten-pmcAfter)/mean(s5new(pmc-w:pmc+w)));
                if( std(sum(blkPmcIP1))/mean(sum(blkPmcIP1))>=0.5)...
                        || (pmcInten-pmcAfter)/mean(s5new(pmc-w:pmc+w))>0.5
                    fprintf('\ncorrect?\n');MCindict(m)=4;
                else fprintf('\nwrong?\n');MCindict(m)=inf;
                end
            end
            
        elseif (pmcEnvlpBefore<pmcEnvlp && pmcEnvlp>pmcEnvlpAfter) &&...
                ~(pmcBefore<pmcInten && pmcInten>pmcAfter)
            fprintf('\n%d a peak in envelop, but not in intensity',pmc);MCindict(m)=1;
        else fprintf ('\n');
        end
        
    end
    
    
else fprintf('MC is empty\n');
end

