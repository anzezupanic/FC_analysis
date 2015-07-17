% function_mmdTestBoot.m
% Maximum mean discrepancy multivariate two sample test
% for testing whether two multivariate distributions are different. It
% gives good results for testing the difference between distribution with
% low number of samples and high dimensionality, but also 
% for distributions with high sample size and lower dimensionality (2-15
% dimensions). the problem with having a large sample size is, that you
% will find significance even for a very small difference between samples, 
% and that (since MMD does not only look at mean, but also variance) it is in practice
% difficult to assign a (biologically) meaningful difference in the distributions.
%
% This is modified code, the original was made available by Gretton et al
% on the following webpage: (downloaded in April 2015)
% http://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm#GreEtAl12
% Author of modifications: Anze Zupanic
%
% The function does pairwise comparison between all distribution in the
% input and computes the MMD test statistics for each. The statistical
% significance is calculated with the bootstrap test.
%
% INPUTs:
% distributions: this is a cell array with each cell representing one
%   distribution (distributions need to be the same size), with samples in
%   rows and variables in columns (for microarrays gene expression values
%   would be in columns, for cytometry the different fluoresences)
% kernelSize: [0,inf]. Gretton et al suggest kernelSize to be equal to the
%   median distance between aggregate samples, but this does not seem to
%   effect calclation much (for large samples). 
% bootSize: the number of bootstrapping runs for calculating significance
% alpha: significanceLevel (normally 0.05)
%
% OUTPUTs:
% solutionTest: in first column you get the MMD values for each comparison,
%   in the second column you get the threshold value for significance and
%   in the third column you get the pvalue
% solutionBoot: in each column all bootstrapped MMD values for asingle
%   comparison are given

function [solutionTest, solutionBoot] = function_mmdTestBoot(distributions, kernelSize, bootSize, alpha)

% a table of all pairs to be compared (all pairs from distributions)
numDist = size(distributions, 2);
pairwise = [];
for ii = 1:numDist-1
    pairwise = [pairwise; ones(numDist-ii,1)*ii, ((ii+1):numDist)'];
end

% prepare storage for output
solutionTest = [];
solutionBoot = [];

%progressBar
progressbar()

for ii = 1:size(pairwise,1)
    i1 = pairwise(ii,1);
    i2 = pairwise(ii,2);
    X = distributions{i1};
    Y = distributions{i2};
    
    % compute the kernel statistics
    m=size(X,1);
    K = rbf_dot(X,X,kernelSize);
    L = rbf_dot(Y,Y,kernelSize);
    KL = rbf_dot(X,Y,kernelSize);
    testStat = 1/m * sum(sum(K + L - KL - KL'));
    
    % bootstrapping
    Kz = [K KL; KL' L];
    MMDarr = zeros(bootSize,1);
    for whichSh=1:bootSize
        [notUsed,indShuff] = sort(rand(2*m,1));
        KzShuff = Kz(indShuff,indShuff);
        K = KzShuff(1:m,1:m);
        L = KzShuff(m+1:2*m,m+1:2*m);
        KL = KzShuff(1:m,m+1:2*m);     
        MMDarr(whichSh) = 1/m * sum(sum(K + L - KL - KL'));
    end
    MMDarr = sort(MMDarr);
    thresh = MMDarr(round((1-alpha)*bootSize));
    solutionTest = [solutionTest; testStat, thresh, sum(testStat<MMDarr)/bootSize];
    solutionBoot = [solutionBoot, MMDarr];
    
    progressbar(ii/size(pairwise,1))
end
    

