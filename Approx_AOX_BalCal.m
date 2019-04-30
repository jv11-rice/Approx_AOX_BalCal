%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize the workspace
clc;
clearvars;
close all;
workspace;
disp('Copyright 2017 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       USER INPUT SECTION
%
out = Approx_AOX_GUI;
if out.cancel == 1
    return
end
%TO SELECT Algebraic Model                         set balCal_FLAG = 1;
%TO SELECT Algebraic and GRBF Model                set balCal_FLAG = 2;
balCal_FLAG = out.grbf;
%
load(out.savePathapp,'-mat');
xapprox = load(out.coeffPath);

%
%
%csvread Read a comma separated value file.
%Adapted from published MathWorks documentation by Mokhtarzadeh, 042516
%
%   M = csvread('FILENAME',R,C,RNG) reads a comma separated value formatted
%   file FILENAME, starting at row R and column C where R=0 and C=0
%   specifies the first value in the file (a zero based index), and reads
%   only the range specified by RNG using spreadsheet notation as in
%   RNG = 'A1..B2' where A1 is the upper-left corner of the data to be
%   read and B2 is the lower-right corner.
%   The file can only contain numeric values.
%
%                       END USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%
% Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        APPROXIMATION       AJM 6/29/17           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
dimFlag = length(excessVecapprox(1,:));
%    load(out.savePathapp,'-mat');
%
% num of data points
numptsapprox = length(excessVecapprox);
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Algebraic Calculations %%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%natural zeros (also called global zeros) are simply set to zero for approximations

globalZerosapprox = mean(natzerosapprox);
%ajm    globalZerosapprox = zeros(length(excessVecapprox(1,:)),1);


%%% make an array out of the globalZerosapprox vector
for i=1:numptsapprox
    globalZerosAllPointsapprox(i,:) = globalZerosapprox;
end
%%%

%%%% 6_14_18 ajm 
matrixcolumnlabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}; 
matrixcolumnlabels(:) = strrep(matrixcolumnlabels(:),'''',''); %get rid of single quotes 

loadlist = {'N1','N2','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'}; 
voltagelist = {'R1','R2','R3','R4','R5','R6','rPLM','rPCM','rMLM','rMCM'}; 
%%%%

%% Subtract the Global Zeros from the Inputs %%%%%%%%%%

for k=1:dimFlag
    
    dainputsapprox(:,k) = excessVecapprox(:,k)-globalZerosAllPointsapprox(:,k);
    
    dalzapprox(:,k) = globalZerosAllPointsapprox(:,k)-globalZerosAllPointsapprox(:,k);
    
end
%%%%%%%%%%%%

% dainputsapprox

% dalzapprox
%numptsapprox


%%
%% Build the Algebraic Model
%%

n(1) = 2*dimFlag*(dimFlag+2);
n(2) = dimFlag*(dimFlag+3)/2;
n(3) = dimFlag;
%model_FLAG = find(n==size(xapprox,1)-1); 
model_FLAG = find(n==size(xapprox,1)); 


%% Full Algebraic Model
if model_FLAG == 1
    nterms = 2*dimFlag*(dimFlag+2);
end

%% Truncated Algebraic Model
if model_FLAG == 2;
    nterms = dimFlag*(dimFlag+3)/2;
end

%% Linear Algebraic Model
if model_FLAG == 3;
    nterms = dimFlag;
end



% Call the Algebraic Subroutine

comGZapprox= zeros(nterms,1);


for i=1:dimFlag
    biggee(:,i) = 0;
end


%[comINapprox,comLZapprox,comGZapprox]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsapprox,0,1,dainputsapprox,dalzapprox,biggee);

%model_FLAG
%nterms
%dimFlag
%numptsapprox
%dainputsapprox
%dalzapprox
%biggee

%series = ones(numptsapprox,1); 

%[comINapprox,comLZapprox]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsapprox,series,1,dainputsapprox,dalzapprox,biggee);

[comINapprox,comLZapprox,comGZapprox]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsapprox,0,0,dainputsapprox,dalzapprox,biggee);
    

    filename = 'AOX_APPROX_comINapprox.csv';
    Z = comINapprox; 
    dlmwrite(filename,Z,'precision','%.16f');       

    filename = 'AOX_APPROX_dainputsapprox.csv';    %AJM 4_24_19
    dlmwrite(filename,dainputsapprox,'precision','%.16f');   


%xapprox

%LOAD APPROXIMATION
%define the approximation for inputs minus global zeros
aprxINapprox = (xapprox'*comINapprox)';        %to find ?? AJM111516

for m=1:length(aprxINapprox)
    ALG_GLOBAL_APPROX(m,:) = aprxINapprox(m,:);
end


%%%%%%
disp(' ');
disp('%%%%%%%%%%%%%%%%%');
disp(' ');
disp('ALG MODEL APPROXIMATION RESULTS: Check ALG_GLOBAL_APPROX.csv in file');

  
    filename = 'ALG_GLOBAL_APPROX.csv';
    dlmwrite(filename,ALG_GLOBAL_APPROX,'precision','%.16f');         
    


%%%%%%%%%%%%%%%%%%%%
%
%%%%    GRBF Section    %%%%%%%
%
%%%%%%%%%%%%%%%%%%%%

%
% Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    RBF SECTION FOR APPROXIMATION     AJM 6/29/17                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if balCal_FLAG == 2
    
    load(out.savePathcal,'-mat');
    GRBF_coeffapprox =  load(out.grbfcoeffPath);
    GRBF_w = load(out.grbfwPath);
    centerIndices =  load(out.grbfcPath);
    
    % GRBF parameters
    numBasis = length(GRBF_coeffapprox);
    ALG_GRBF_GLOBAL_APPROX = ALG_GLOBAL_APPROX;
    
    %natural zeros (also called global zeros) are simply set to zero for approximations
    
    globalZerosCalib = mean(natzerosCalib);
    
    %%% make an array out of the globalZerosCalib vector
    
    %    for i=1:numptsCalib
    for i=1:length(excessVecCalib)
        globalZerosAllPointsCalib(i,:) = globalZerosCalib;
    end
    %%%


   
for k=1:dimFlag
    dainputscalib(:,k) = excessVecCalib(:,k)-globalZerosAllPointsCalib(:,k);
end    
    
    for u=1:numBasis
        
        for s=1:dimFlag
            
            
            centerIndexLoop(s) = centerIndices(u,s);
%            w(s) = GRBF_w(u,s);
            
            
            for r=1:length(excessVecapprox(:,1))
                etaapprox(r,s) = dot(dainputsapprox(r,:)-dainputscalib(centerIndices(u,s),:),dainputsapprox(r,:)-dainputscalib(centerIndices(u,s),:));
            end
            
            
            rbfINminGZapprox(:,s)=exp(etaapprox(:,s)*log(abs(GRBF_w(u,s))));
rbfc_INminGZapprox(:,s) = GRBF_coeffapprox(u,s)*rbfINminGZapprox(:,s); % ajm 6_7_18
%            rbfc_INminGZapprox(:,s) = GRBF_coeffapprox(s)*rbfINminGZapprox(:,s);
            
        end
        
        
        
        %UPDATE THE RESIDUAL
        
        %update the approximation
        
        ALG_GRBF_GLOBAL_APPROX = ALG_GRBF_GLOBAL_APPROX  + rbfc_INminGZapprox;
        
        rbfc_INminGZ_Histapprox{u} = rbfc_INminGZapprox;  % temp ajm 6_7_18
   
   rbf_etaapprox_Hist{u} = etaapprox;   % temp ajm 6_7_18
   
            rbf_excessveccalib_center_Hist{u} = dainputscalib(centerIndices(u,:),:);  % temp ajm 6_7_18
    end
    

total_rbfc_INminGZ_approx = rbfc_INminGZ_Histapprox{1};            % temp ajm 6_7_18
    
    for u=2:numBasis
total_rbfc_INminGZ_approx = total_rbfc_INminGZ_approx + rbfc_INminGZ_Histapprox{u};            % temp ajm 6_7_18
    end
    
%%%%%%
disp(' ');
disp('%%%%%%%%%%%%%%%%%');
disp(' ');
disp('ALG + GRBF MODEL APPROXIMATION RESULTS: Check ALG+GRBF_GLOBAL_APPROX.csv in file');


    filename = 'ALG+GRBF_GLOBAL_APPROX.csv';
    dlmwrite(filename,ALG_GRBF_GLOBAL_APPROX,'precision','%.16f');         
  
%%%%%%%%
    

    %
    %End Approximation Option
end


disp('  ')
disp('Calculations Complete.')

% Tidy up the Workspace
clearvars -except ALG_GLOBAL_APPROX ALG_GRBF_GLOBAL_APPROX

%
% Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
%

%toc

