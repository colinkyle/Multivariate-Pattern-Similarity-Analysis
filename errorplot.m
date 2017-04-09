function [varargout] = errorplot(Data,varargin)
%Bar plot with errorbars and shows significance for tests between bars
%
%[varargout] = errorplot(Data,'xvalues',[1:size(Data,2)],'verbose',false,
%       'color',[0 .55 .65],'testtype','paired','alpha',.05,)
%
%   Applies bonferroni correction for multiple comparisons
%
%   Example:
%   Data = randn(10,10)+[zeros(10,5),1*ones(10,2),2*ones(10,3)]
%   errorplot(Data,'testtype','unpaired','color',[.5,.5,.5],'verbose',true)
parser = inputParser;
parser.addRequired('Data',@ismatrix);
parser.addOptional('verbose',false,@islogical);
parser.addOptional('color',[0 .55 .65],@ismatrix);
parser.addOptional('xvalues',[1:size(Data,2)],@ismatrix);
parser.addOptional('testtype','paired',@isstr);
parser.addOptional('alpha',0.05,@isfloat);
parser.parse(Data,varargin{:})

alpha = parser.Results.alpha;
color = parser.Results.color;
testtype = parser.Results.testtype;
verbose = parser.Results.verbose;
xvalues = parser.Results.xvalues;

Dm = nanmean(Data);
Dse = nanstd(Data)./sqrt(size(Data,1));
hand = bar(xvalues,Dm,'FaceColor',color,'LineWidth',3);
hold on

drawerrorbars(xvalues,Dm,Dse)

height0 = 1.1*max([Dm+Dse]);
dh = 0.07*max([Dm+Dse]);
nComparisons = size(Data,2);
[h,p] = ttest(Data,0,'alpha',alpha/nComparisons);
p = p*nComparisons;

if verbose
    for ii = 1:length(h)
        if h(ii)==1
            % put dots here
            if p(ii)<=.01
                plot([xvalues(ii)-.1,xvalues(ii)+.1],[0,0],'k*','linewidth',1)
            elseif p(ii)<=.05
                plot(xvalues(ii),0,'k*','linewidth',1)
            end
            text(xvalues(ii),.5*dh,num2str(p(ii),'%2.3f'),'HorizontalAlignment','center','fontsize',12,'fontweight','bold')
        end
    end
else
    for ii = 1:length(h)
        if h(ii)==1
            % put dots here
            if p(ii)<=.01
                plot([xvalues(ii)-.1,xvalues(ii)+.1],[0,0],'k*','linewidth',1)
            elseif p(ii)<=.05
                plot(xvalues(ii),0,'k*','linewidth',1)
            end
        end
    end
end

nComparisons = (size(Data,2)*(size(Data,2)-1))/2;
for i1 = 1:(size(Data,2)-1)
    for i2 = (i1+1):size(Data,2)
        if strcmp(testtype,'paired')
            [h,p] = ttest(Data(:,i1),Data(:,i2),'alpha',alpha/nComparisons);
        else
            [h,p] = ttest2(Data(:,i1),Data(:,i2),'alpha',alpha/nComparisons);
        end
        %%
        p = p*nComparisons;
        if verbose
            if h==1
                plot([xvalues(i1),xvalues(i2)],[height0,height0],'k','linewidth',1)
                plot([xvalues(i1),xvalues(i1)],[height0,height0-.3*dh],'k','linewidth',1)
                plot([xvalues(i2),xvalues(i2)],[height0,height0-.3*dh],'k','linewidth',1)
                if p<=.01
                    plot(mean([xvalues(i1),xvalues(i2)])-.1,height0+.25*dh,'k*','linewidth',1)
                    plot(mean([xvalues(i1),xvalues(i2)])+.1,height0+.25*dh,'k*','linewidth',1)
                elseif p<=.05
                    plot(mean([xvalues(i1),xvalues(i2)]),height0+.25*dh,'k*','linewidth',1)
                end
                height0 = height0+dh;
                text(mean([xvalues(i1),xvalues(i2)]),height0-1.5*dh,num2str(p,'%2.3f'),'HorizontalAlignment','center','fontsize',10,'fontweight','bold')
            end
        else
            if h==1
                
                plot([xvalues(i1),xvalues(i2)],[height0,height0],'k','linewidth',1)
                plot([xvalues(i1),xvalues(i1)],[height0,height0-.3*dh],'k','linewidth',1)
                plot([xvalues(i2),xvalues(i2)],[height0,height0-.3*dh],'k','linewidth',1)
                if p<=.01
                    plot(mean([xvalues(i1),xvalues(i2)])-.1,height0+.25*dh,'k*','linewidth',1)
                    plot(mean([xvalues(i1),xvalues(i2)])+.1,height0+.25*dh,'k*','linewidth',1)
                elseif p<=.05
                    plot(mean([xvalues(i1),xvalues(i2)]),height0+.25*dh,'k*','linewidth',1)
                end
                height0 = height0+dh;
            end
        end
        
    end
end
varargout{1}=height0;
varargout{2}=dh;
varargout{3}=hand;
xlim([.3 size(Data,2)+.8])
box off
end

function drawerrorbars(xvalues,Dm,Dse)
barwid = .2*.8*(mean(diff(xvalues))./2);
for i = 1:numel(xvalues)
    plot([xvalues(i),xvalues(i)],[Dm(i)-Dse(i)/2,Dm(i)+Dse(i)/2],'k','linewidth',3)
    plot([xvalues(i)-barwid,xvalues(i)+barwid],[Dm(i)-Dse(i)/2,Dm(i)-Dse(i)/2],'k','linewidth',3)
    plot([xvalues(i)-barwid,xvalues(i)+barwid],[Dm(i)+Dse(i)/2,Dm(i)+Dse(i)/2],'k','linewidth',3)
end
end

