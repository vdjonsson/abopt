clear all 
clc
close all 

%% This is the code for analyzing the NIH data for saturated data.  

basepath = '/Users/vjonsson/biodata/hivabcomb/';
fname = strcat(basepath, 'ablistone.csv');

ablist = importdata(fname,'\n',200)
abname = ablist{1}; 

filesloaded = cell(1);

fpath = '/Users/vjonsson/biodata/hivabcomb/';
mp = 'matrix/';

matpath = strcat(fpath, mp);

cd(matpath) 

files = dir('*.csv');
jjn =1; 
jjs =1;

strneut = strcat('MatrixNeut_', abname,'.csv');
strseq = strcat('MatrixSeq_', abname,'.csv');

% Import the data 
for i=1:length(files)
    eval(['load ' files(i).name ' -ascii']);
    add =0;
    fname  = files(i).name;
    
    if  strfind(fname, strneut)
        B= importdata(fname); 
        filesloaded{jjn} = fname;
        jjn = jjn+1;
    end
    if  strfind(fname, strseq)
        A= importdata(fname);
        jjs = jjs+1;
    end    
end


% File output 
fp = '/Users/vjonsson/biodata/hivabcomb/'
cvpath = 'cv/';

solfp ='fit/';
figsp = 'fig/';
fext ='.csv';

sizealldata = size(A,1);

indexp = 1:1:sizealldata; 

% First put saturation points in a different matrix as well saturation level points  
indmax = min(find(B==max(B)))-1
satindex = indmax;

AT = A(1:satindex,:);
BT = B(1:satindex);

ATS = A(satindex+1:sizealldata,:);
BTS = B(satindex+1:sizealldata);

lam_1 = 10:10:100;
lam_2 = 10:10:60;

numpartitions = 5; 

ii = 1;

for i=1:length(lam_2)
  for j=1:length(lam_1)
    lam2 = lam_2(i);
    lam1 = lam_1(j);

    if (lam1 + lam2) <= 100
        lam3 = 100 - lam1 -lam2;

        if lam2 >= 0 && lam1 >=0
            lenl2 = length(lam_2);
            lamarrsat(ii,:) = [lam1, lam2, lam3];  
            ii = ii+1;

        end
    end
  end
end


lamarrsat

weights = 0; % fix this 

% First do the saturated lasso and cross validation generating a good
% lambda
[lamparamsat, toterrsat,figcv] = satlassocrossval(A, B, satindex, abname, numpartitions, lamarrsat,weights)

% xfn = strcat(fp, solfp, 'LAMBDA_', abname, fext);
% csvwrite(xfn, lamparamsat);

% 
% % toterrsat is the lowest total error: model plus estimation error picked from 
% % % cross validation 
% allblue = 1; 
% [x,ysat] = augmentedlassosolution(AT, BT, ATS, BTS, lamparamsat);
% fig = graphlassosolution(A, B, A, indexp, x, ysat, lamparamsat, abname, 'SAT LASSO,', allblue);
% totalerrsatmodel = erroraugmentedlasso(AT, BT, ATS, BTS, ysat) + modelerrorsatlasso(AT,BT,ATS,BTS, ysat)
% 
% % write to file 
% xfn = strcat(fp, solfp, 'XS_', abname, fext);
% csvwrite(xfn, ysat);
% 
% fitfig = strcat(fp,figsp,'FS_', abname,'.jpg');
% saveas(fig, fitfig);
% 
% figcvname = strcat(fp,figsp,'FCV_', abname,'.jpg');
% saveas(figcv, figcvname);


% % % First do the regular lasso 
% [lamparam, toterr] = lassosolutioncrossval(A, B, satindex, abname, numpartitions, lamarr);
%%
[x,y] = satlassosolution(AT, BT, ATS, BTS, lamparamsat,weights);
fig = graphlassosolution(A, B, A, indexp, x, y, lamparamsat, abname, 'LASSO,',1);

% esterr = erroraugmentedlasso(AT, BT, ATS, BTS, y);
% totalerrmodel = erroraugmentedlasso(AT, BT, ATS, BTS, y) + modelerrorsatlasso(AT,BT,ATS,BTS, y)

% % write to file 
xfn = strcat(fpath,solfp, 'XW_', abname, fext);
csvwrite(xfn, y);
% 
fitfig = strcat(fpath,figsp,'FL_', abname,'.png');
saveas(fig, fitfig);

beep
beep 
beep 

