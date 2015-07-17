% viSNE_MMD.m
% this script calculates the maximum mean discrepancy between loaded viSNE
% submaps. Depending on the size of the datasets the algorithm runtime is 
% ca. 1-60 minutes. We recommend using at most 5000 cells/particles per
% sample to keep runtime reasonbale low. All samples must be of the sme
% size, otherwise comparison between MMD values is pointless. 

% load data 
%   field samples:          environment_biological_similarity_3000_SNE.mat
%   lab samples:            lab_timesandtreatment_3000_SNE.mat
%   single species samples: 
load('lab_timesandtreatment_3000_SNE.mat')

% load the viSNE submaps  - for the lab samples, the viSNE variables can be
% found in the 54 and 55 columns of sessionData. For other datasets, check
% the third column gated to find the bh-SNE variable.
for ii = 1: size(gates,1)
    index = gates{ii,2};
    distributions{ii}=sessionData(index,[54,55]);
end

% Check the format of distributions (cell array). The MMD function needs
% this format to run (The file funtion_mmdtestBoot.m needs to be in the
% same folder as this script). 
% In this call we use the following values:
%    kernelSize: 1 Gretton et al suggest kernelSize to be equal to the
%       median distance between aggregate samples, which in our case is approximately one, 
%       but zhis choice does not affect calculations much when you have a
%       high enough number fo samples
%   bootSize: 1 if you want to use bootstrapping to calculate the
%       significance of the MMD statistics (i.e., to calculate whether two
%       samples are significantly different) use at least 100. 
%   alpha: 0.05 - significance level 
[solutionTest, solutionBoot] = function_mmdTestBoot(distributions, 1, 1, 0.05);

% in MMDvalues you find the maximum discrepancy for each pairwise
% comparison performed
MMDvalues = solutionTest(:,1);
