close all; 
clear lamarr;
clear lamparam
clear lam_1;
clear lam_2; 
clear lam_3;
clear lam1;
clear lam2; 
clear lam3;
clc; 

% File output 
fp = '/Users/vjonsson/data/eclipse/workspace/statshiv/fileparse/statsab/';
cvpath = 'cv/';

solfp ='fits/';
figsp = 'figs/';
fext ='.csv';

abname = '3BC176';
A = MatrixSeq3BC176;
B = MatrixNeut3BC176;

sizealldata = size(A,1);

indexp = 1:1:sizealldata; 

% First put saturation points in a different matrix as well saturation level points  
indmax = min(find(B==max(B)))-1;
satindex = indmax;

AT = A(1:satindex,:);
BT = B(1:satindex);

ATS = A(satindex+1:sizealldata,:);
BTS = B(satindex+1:sizealldata);

lam_1 = 0:10:100;
lam_3 = 0:10:100; 

numpartitions = 5; 

% lamparam = lassosolutioncrossval(fullA, fullB, satindex, abname, numparts, lamarr)
ii = 1;

for i=1:length(lam_3)
  for j=1:length(lam_1)
    lam3 = lam_3(i);
    lam1 = lam_1(j);
    
    if (lam1 + lam3) <= 100
        lam2 = 100 - lam1 -lam3;

        if lam2 > 0 && lam1 >0 && lam3 > 0  

            lenl2 = length(lam_3);

            lamarr(ii,:) = [lam1, lam2, lam3];  
            ii = ii+1;

            lamparam = augmentedlassosolutioncrossval(A, B, satindex, abname, numpartitions, lamarr)
        end
    end
  end
end

[x,y] = augmentedlassosolution(AT, BT, ATS, BTS, lamarr);

fig = graphlassosolution(A, B, A, indexp, x, y, lamparam, abname)
 
mse = erroraugmentedlasso(AT, BT, ATS, BTS, y)

% write to file 
xfn = strcat(fp, solfp, 'X_', abname, fext);
csvwrite(xfn, y);

% write to file 
fitfig = strcat(fp,figsp,'FIG_', abname,'.jpg') 
saveas(fig, fitfig);

        