function [ uSq ] = dir_watsonusq( varargin )
% Calculates the Watson Usq statistic for one or two sample test (values in deg)
% Watson Usq is a non-parametric test for verifying divergence between two circular distributions or
% for the divergence of a single circular distribution from a perfect circle.
%
%       uSq = dir_watsonusq(dir);
%       uSq = dir_watsonusq(dir1,dir2);
%       
% Test either one directional rate map against a cirular distribution, or
% two rate maps against each other.

% NB have checked the 2 sample code against a worked example - it's correct and the same for the one
% sample code it matches the results in Zar
%
% Code by Caswell Barry, 2008.


% -------------------------------------------------------------------------------------------------
% --- PARSE ARGUMENTS -----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
switch nargin
    case 1
        %Do a single sample test - distribution vs a circle

        ui=sort(varargin{1}(:)./360);

        n=length(ui);
        sumUi=sum(ui);
        sumUiSq=sum(ui.^2);
        meanUi=mean(ui);
        sumiUi=sum((1:n)'.*ui);

        uSq=sumUiSq - sumUi.^2/n - (2/n)*sumiUi + (n+1)*meanUi + n/12;


    case 2
        %Do a two sample test - distribution vs each other
        sample1=sort(varargin{1}(:));
        sample2=sort(varargin{2}(:));
        clear varargin;

        %Need to get unique numbers and rank the numbers then normalise the ranks
        rank1=1:length(sample1);
        rank2=1:length(sample2);

        [unique1, ind1]=unique(sample1);
        [unique2, ind2]=unique(sample2);
        normRank1=rank1(ind1)./length(sample1);
        normRank2=rank2(ind2)./length(sample2);


        %Then need to compare the two lists of numbers - to get to this need several steps
        %First find unique list for both samples and also record the number of times each number
        %appears
        [uniqueAll, nInstancesAll]= ilUniqueInst([sample1;sample2]);

        finalList1=ones(length(uniqueAll),1)*nan;
        finalList1=[0;finalList1]; %Add a 0 to the start - remember to take off later
        finalList2=finalList1;

        for n=1:length(uniqueAll)
            if any(unique1==uniqueAll(n))
                finalList1(n+1)=normRank1(unique1==uniqueAll(n));
            else
                finalList1(n+1)=finalList1(n);
            end

            if any(unique2==uniqueAll(n))
                finalList2(n+1)=normRank2(unique2==uniqueAll(n));
            else
                finalList2(n+1)=finalList2(n);
            end
        end

        %Find difference and square and remove leading zero
        d=(finalList1-finalList2);
        d=d(2:end);
        dSq=d.^2;


        %Final few calcs
        n1=length(sample1);
        n2=length(sample2);

        %Finally calc Usq
        uSq=((n1*n2)/(n1+n2).^2) * (sum(nInstancesAll.*dSq) - (sum(nInstancesAll.*d).^2)/(n1+n2));


    otherwise
        %Test doesn't make sense for other number of samples
        error('Watson test can only be run for one or two samples');
end



function [uniqueList, nInstances]= ilUniqueInst(list)
% From a numerical list finds the unique numbers (which built in function does anyway) but also
% returns a second vector indicating how many times each unique number was in the original list.
% Uses one of the functions of unique.m to do this. Given B=unique(A) - unique.m will also return a
% vector i such that B=A(i) if a number is represented multiple times in A then i will
% always index the last instance of it.
list=sort(list);
[uniqueList, ind2Org]=unique(list);
nInstances=ind2Org-[0;ind2Org(1:end-1)];
