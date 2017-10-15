function [h, p, table, ptable] = chi_squared(data, varargin)
% CHI_SQUARED performs chi-squared test.
%
% [h, p, table, ptable] = chi_squared({x1, x2}, ...) takes categorical
% variables, vectors, string arrays, or cell arrays of strings.
%
% If x1 is a numeric vector taking all values 1,2,...,M, and x2 is a
% numeric vector taking all values 1,2,...N, then table is an M-by-N
% matrix with table(i,j) equal to the number of elements with x1=i and x2=j.
%
% All inputs must be the same length.
%
% The function then conducts a chi-squared test on the table, and returns the
% test p-value in p.
%
% It further conducts post-hoc test by getting the multiple comparison corected p-value
% of the adjusted residuals in each cell of the table. This is based on 
%
%   Beasley, T. Mark, and Randall E. Schumacker. "Multiple regression approach to 
%   analyzing contingency tables: Post hoc and planned comparison procedures." 
%   The Journal of Experimental Education 64.1 (1995): 79-93.
%
% The ptable is a logical array indicating whether the adjusted p-value of the
% adjusted residual is significant.
%
% Inputs:
%   - data: cell of input vectors, should all be the same length.
% 
% varargin:
%   - 'fadj'  : Gives the name of multiple comparison procedure. Valid ones includes
%              'bonferroni' (Default), and 'fdr'.
%   - 'alpha' : p-value threshold for significance, for both initial chi2 test and the
%               post-hoc tests. Default is 0.05.
%
% Output:
%   - h:     1 if ch2 shows the categorical variables are not independent. 0 otherwise.
%   - p:     pvalue of the initial chi2 test.
%   - table: Contingency table from the supplied data vectors
%   - ptable: Table indicating whether each cell is significantly different than expected, corresponds
%             to table. Values can be 0 for not different than expected, 1 means more frequent than
%             expected, -1 means less frequent than expected.
%
% See also CROSSTAB
    
    funcInputs = parseMyInputs(data, varargin{:});

    % First make the contingency table -- similar to CROSSTAB
    sz = zeros(1, numel(data));
    nonan = [];
    M = [];
    for j = 1:numel(data)
        % g1=which group each element belongs to
        % g2=value of the different groups in this category
        [g1, g2] = grp2idx(data{j}); 
        ng = size(g2, 1);   % number of groups for each category
        sz(j) = ng;
        n = length(g1);     % number of elements in this category vector
        if (j==1)
            n1 = n;         % n1=number of total elements every vector should have
            nonan = ~isnan(g1);
            M = zeros(n, numel(data));  % Columns of M is the groupIdx for each category
            M(:,1) = g1;
        elseif (n ~= n1)
            error(message('chi-squared: InputSizeMismatch'));
        else
            nonan = nonan & ~isnan(g1); % keep track what elements have labels in all categories
            M(:,j) = g1;
        end
    end

    M = M(nonan, :);    % only keep elements with labels in all categories
    n = size(M, 1);     % total elements

    table = zeros(sz);
    % make table row element by element.
    for k = 1:n
        N = num2cell(M(k,:));
        table(N{:}) = table(N{:}) + 1;
    end

    % remove degenerate dimensions, i.e. categories with only 1 group
    if any(sz==1)
        sz = sz(sz>1);
        T = reshape(table, sz);
    else
        T = table;
    end

    % chi-squared test procedures
    % crazy-way to get the expected value..
    expected = zeros(sz);   
    expected(:) = n;
    szv = sz;
    permv = [(2:length(sz)), 1];
    for j = 1:numel(data)
        sz1 = szv(1);
        sz2 = prod(szv(2:end));
        frac = sum(reshape(T, sz1, sz2), 2)/n;
        frac = reshape(repmat(frac,1,sz2),szv);
        expected = expected.*frac;
        expected = shiftdim(expected, 1);
        szv = szv(permv);
        T = shiftdim(T, 1);
    end

    % chi-squared value
    chi2 = (T - expected).^2./expected;
    chi2 = sum(chi2(:));
    df = prod(sz) - (1+sum(sz-1));
    p = gammainc(chi2/2, df/2, 'upper');    % see stat toolbox chi2pval
    h = (p<=funcInputs.alpha);              % chi2 omnibus -- reject null hypothesis?

    % now do residual testing
    ptable = (T - expected)./sqrt(expected);    % standardized residuals, or z-scores
    ptable = 2.*(1-normcdf(abs(ptable))).*sign(ptable); % magnitude is two-tailed pvalue.
    if strcmp(funcInputs.fadj, 'bonferroni')
        alpha_adj = funcInputs.alpha/numel(T);  % adjusted alpha-value
        ptable = (abs(ptable)<=alpha_adj).*sign(ptable);
    elseif strcmp(funcInputs.fadj, 'fdr')
        ptable(:) = mafdr(abs(ptable(:)), 'bhfdr', true).*sign(ptable);
        ptable = (abs(ptable)<=funcInputs.alpha).*sign(ptable);
    end
end

function funcInputs = parseMyInputs(data, varargin)
    funcInputs = inputParser;
    addRequired(funcInputs, 'data', @(x) numel(x)>=2);  % need at least two factors
    addParameter(funcInputs, 'fadj', 'bonferroni', @(x)strcmp(x,'bonferroni') | strcmp(x,'fdr'));
    addParameter(funcInputs, 'alpha', 0.05, @(x) isscalar(x) & x>0);
    parse(funcInputs, data, varargin{:});
    funcInputs = funcInputs.Results;
end

