function STATS=testt(x1,x2,varargin)
% Student's t test for unpaired or paired samples.
% This file is applicable for equal or unequal sample sizes; for paired or 
% unpaired samples. When the test is unpaired, the Fisher-Snedecor F-test is
% performed to assess the equality of variance. If variances are not equal,
% Satterthwaite's approximate t test is performed.
% 
% Syntax: 	TESTT(X1,X2,TST,ALPHA,TAIL)
%      
%     Inputs:
%           X1 and X2 - data vectors (mandatory). 
%           TST - unpaired (0) or paired (1) test (default = 0).
%           ALPHA - significance level (default = 0.05).
%           TAIL - 1-tailed test (1) or 2-tailed test (2). (default = 1).
%     Outputs:
%           - Fisher-Snedecor test (for unpaired test)
%           - t value.
%           - degrees of freedom.
%           - Confidence interval of means difference (for paired test)
%           - p-value
%           - Power
% 
%      Example: 
% 
%           X1=[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 ...
%           85 86 86 87 87];
% 
%           X2=[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 ...
%           88 88 89 90 90];
% 
%           Calling on Matlab the function: testt
% 
%           Answer is:
% 
% FISHER-SNEDECOR F-TEST FOR EQUALITY OF VARIANCES
%  
%       F       DF_numerator    DF_denominator    p_value
%     ______    ____________    ______________    _______
% 
%     1.5379    24              24                0.29861
% 
% --------------------------------------------------------------------------------
% Variances are equal
% --------------------------------------------------------------------------------
%  
% STUDENT'S T-TEST FOR UNPAIRED SAMPLES
%  
%       t       DF    tail    alpha     p_value      Power 
%     ______    __    ____    _____    _________    _______
% 
%     5.2411    48    1       0.05     1.765e-06    0.99958
% 
% STATS=TESTT(...) returns a structure with all test(s) statistics
% 
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
% 
% To cite this file, this would be an appropriate format:
% Cardillo G. (2006). Student t-Test for unpaired or paired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/12699

%Input error handling
p = inputParser;
addRequired(p,'x1',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
addRequired(p,'x2',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
addOptional(p,'tst',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'tail',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==1 || x==2));
parse(p,x1,x2,varargin{:});
tst=p.Results.tst; alpha=p.Results.alpha; tail=p.Results.tail;
clear p
if tst==1
    assert(length(x1)==length(x2),'Warning: x1 and x2 must have the same length')
end
tr=repmat('-',1,80);

switch tst
    case 0 %unpaired test
        n=[length(x1) length(x2)]; %samples sizes
        m=[mean(x1) mean(x2)]; %samples means
        v=[var(x1) var(x2)]; %samples variances
        %Fisher-Snedecor F-test
        if v(2)>v(1) 
            v=fliplr(v);
            m=fliplr(m);
            n=fliplr(n);
        end
        F=v(1)/v(2); %variances ratio
        DF=n-1;
        p = fcdf(F,DF(1),DF(2)); %p-value
        p = 2*min(p,1-p);
        if nargout
            STATS.Fvalue=F;
            STATS.DFn=DF(1);
            STATS.DFd=DF(2);
            STATS.FPvalue=p;
        end
        %display results
        disp('FISHER-SNEDECOR F-TEST FOR EQUALITY OF VARIANCES')
        disp(' ')
        disp(array2table([F DF(1) DF(2) p],'VariableNames',{'F' 'DF_numerator' 'DF_denominator' 'p_value'}));
        disp(tr)
        if p<alpha %unequal variances (Behrens-Welch problem)
            fprintf('Variances are different: Behrens-Welch problem\n')
            disp(tr)
            disp(' ')
            %Satterthwaite's approximate t test
            a=v./n; b=sum(a);
            denom=sqrt(b);
            gl=b^2/sum(a.^2./(n-1));
            disp('SATTERTHWAITE''S APPROXIMATE T-TEST FOR UNPAIRED SAMPLES')
            disp(' ')
        else %equal variances
            fprintf('Variances are equal\n')
            disp(tr)
            disp(' ')
            gl=sum(n)-2; %degrees of freedom
            s=sum((n-1).*v)/(sum(n)-2); %combined variance
            denom=sqrt(sum(s./n));
            disp('STUDENT''S T-TEST FOR UNPAIRED SAMPLES')
            disp(' ')
        end
        dm=diff(m); %Difference of means
        clear H n m v a b s %clear unnecessary variables
    case 1 %paired test
        disp('STUDENT''S T-TEST FOR PAIRED SAMPLES')
        disp(' ')
        disp(tr)
        n=length(x1); %samples size
        gl=n-1; %degrees of freedom
        d=x1-x2; %samples difference
        dm=mean(d); %mean of differences
        vc=tinv(1-alpha/tail,gl); %critical value
        ic=[abs(dm)-vc abs(dm)+vc]; %Confidence interval
        denom=sqrt((sum((d-dm).^2))/(n*(n-1))); %standard error of difference
        TBL=table(abs(dm),ic);
        TBL.Properties.VariableNames = {'Mean_of_differences' 'Confidence_interval'};
        disp(TBL)
        clear n d %clear unnecessary variables
end
t=abs(dm)/denom; %t value
p=(1-tcdf(t,gl))*tail; %t-value associated p-value

switch tail
    case 1
        if p>=alpha
            tp = tinv(1-alpha,gl) - t;  %Power estimation.
            Power=1-tcdf(tp,gl);
        else
            tp = t - tinv(1-alpha,gl);  %Power estimation.
            Power=tcdf(tp,gl);
        end
    case 2
        if p>=alpha
            tp1 = tinv(1-alpha/2,gl) - t;  %Power estimation.
            tp2 = t + tinv(1-alpha/2,gl);
            Power=2-tcdf(tp1,gl)-tcdf(tp2,gl);
        else
            tb1=t - tinv(1-alpha/2,gl);  %Power estimation.
            tb2=t + tinv(1-alpha/2,gl);
            Power=1 - (tcdf(tb2,gl)-tcdf(tb1,gl));
        end
end
disp(array2table([t gl tail alpha p Power],'VariableNames',{'t' 'DF' 'tail' 'alpha' 'p_value' 'Power'}))

if nargout
    STATS.tvalue=t;
    STATS.tdf=gl;
    STATS.ttail=tail;
    STATS.tpvalue=p;
    STATS.power=Power;
end

end
