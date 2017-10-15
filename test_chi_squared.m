function tests = test_chi_squared()
    tests = functiontests(localfunctions);
end

% Test data adopted from Table 1 of:
%
%   Beasley, T. Mark, and Randall E. Schumacker. "Multiple regression approach to 
%   analyzing contingency tables: Post hoc and planned comparison procedures." 
%   The Journal of Experimental Education 64.1 (1995): 79-93.
%
%                           Prep
%       BOTH        COGN       CULT     NONE
% EE    59          35          8        48     | 150
% ME    23          44          2        25     | 94
% MM    35          27          20       24     | 106
%       ----------------------------------------|
%       117        106          30       97      350
%
% chi2=37.026, p<0.00001

function test1(testCase)
% PREP values are 0, 1, 2, 3 --> BOTH, COGN, CULT, NONE
% ED values are   0, 1, 2    --> EE, ME, MM
    
    PREP = [zeros(117,1); ones(106, 1); 2*ones(30, 1); 3*ones(97,1)];
    ED =   [zeros(59,1); ones(23, 1); 2*ones(35,1);
            zeros(35,1); ones(44, 1); 2*ones(27,1);
            zeros(8,1);  ones(2,1);   2*ones(20,1);
            zeros(48,1); ones(25,1);  2*ones(24,1)];


    [h, p, table, ptable] = chi_squared({ED,PREP}, 'alpha', 0.05);

    assert(p<0.00001);
    assert(h==1);
    assert(all(all(table==[59, 35, 8, 48; 23, 44, 2, 25; 35, 27, 20, 24])));
    assert(all(all(ptable==[0, 0, 0, 0; 0, 1, 0, 0,; 0, 0, 1, 0])));
end


