function [lamparam, totalerr] = satlassocrossval(fullA, fullB, satindex, abname, numpartitions, lamarr, regtype)

% This function does cross validation on the full data set, including
% saturation points and uses the augmented lasso function that penalizes 
% the saturated data 
% Given number of partitions, partition into training set and validation 
% set.  Each element from these sets are picked randomly, with only 20% of
% the training and validation sets picked ramdomly from the saturated set 


    fnpath = strcat('/Users/vjonsson/Google Drive/biodata/hivstats/matrix/', abname, '/');

    xfn = strcat(fnpath, 'ALLLAM_', abname,'_', regtype, '.csv');
    csvwrite(xfn, lamarr);
    
    A = fullA; 
    B = fullB; 
    bstr = abname; 
    lenlam = size(lamarr,1);
    
    % First separate all the data into good data and saturated data     
    sizealldata = size(fullB,1);
    n = size(fullA,2); %col space, max size of regressors

    for jj=1:lenlam
        
        lamarr(jj,:)
        
        jj;
        lam1 = lamarr(jj, 1);
        lam2 = lamarr(jj, 2);
        lam3 = lamarr(jj, 3);
    
        % Indices for good and saturated data 
        BTi = 1:satindex;
        BTSi = satindex+1:sizealldata;
        
        satval = fullB(satindex+1);
       
        pctrain = 1 -1/numpartitions; 
     
        sizesatdata = length(BTSi);
        sizedata = sizealldata-sizesatdata; 
        
        trainsize = floor(sizealldata*pctrain); 
        valsize = floor(sizealldata - trainsize);
        
        % For each training partition pick off x% points from good data 
        % and y% points from the saturated data, where x and y are
        % calculated depending on the ratio of good and bad data 
        
        pcsatincl = sizesatdata/sizealldata; 
         
        trainsatsize = floor(trainsize*pcsatincl); 
        traingoodsize = trainsize - trainsatsize;
        
        if traingoodsize > length(BTi)
            traingoodsize = traingoodsize -1;
        end
     
        for j=1:numpartitions
            
            trgoodinds = datasample(BTi,traingoodsize,'Replace', false); 
            trsatinds = datasample(BTSi,trainsatsize);  
            
            % Now get validation indices different from training ones
            i= 1; 
            
            trallinds = [trgoodinds, trsatinds];
            
            for k=1:sizealldata
                if isempty(find(trallinds == k))
                    valindexset(i)=k;
                    i= i+1;
                end
            end
            
            % If the valindexset is empty, then just pick random validation
            % set 
            if i > 1
                % From here, pick the random validation index set, it should not matter 
                valinds = datasample(valindexset,valsize);
            else 
                valinds  = datasample(trallinds, valsize);
            end
            
            % Now construct all sets 
            
            %Good training data 
            for i=1:length(trgoodinds)
                AT(i,:) = fullA(trgoodinds(i),:);
                BT(i,1) = fullB(trgoodinds(i));
            end
            
            %Saturated training data 
            for i=1:length(trsatinds)
                ATS(i,:) = fullA(trsatinds(i),:);
                BTS(i,1) = fullB(trsatinds(i));
            end
       
            l = 1;
            ll = 1;
            
            AVS = [];
            BVS =[];
            %Validation data 
            for i=1:length(valinds)
                
                AVV(i,:) = fullA(valinds(i),:);
                BVV(i,1) = fullB(valinds(i));
                
                if valinds(i) < satindex
                    AV(l,:) = fullA(valinds(i),:);
                    BV(l,1) = fullB(valinds(i));
                    valindgood(l) = valinds(i);
                    l = l +1;
                else
                    AVS(ll,:) = fullA(valinds(i),:);
                    BVS(ll,1) = fullB(valinds(i));
                    valindsat(ll) = valinds(i);
                    ll = ll +1;
                end
            end
            
            if isempty(BVS)
                BVS(1) = satval;
            end
            
            
            % concatenate good and saturated data to graph later
            ATT = [AT; ATS];
            
            xt = trallinds;
            xv = valinds;

            msize = 1/length(BT); % size of the good data 

            [y] = satlasso(AT, BT, ATS, BTS, lamarr(jj,:), regtype);

            esterr= erroraugmentedlasso(AV, BV, AVS, BVS, y)
            modelerr = 0;
            
            totalerr = esterr + modelerr
            
            MSE(jj, j) = esterr;
            ME(jj, j) = modelerr; 
            
            Av = AVV*y;
            Ay = ATT*y;
             
         
        end   
         
        CVE(jj) = 1/numpartitions*sum(MSE(jj,:));
        MDE(jj) = 1/numpartitions*sum(ME(jj,:));
        TOTE(jj) = MDE(jj) + CVE(jj);
        
      
    end

    % Pick the y that has the lowest TOTE 
    CVE 
    MDE     
    TOTE
    minTOTE = min(TOTE);   
    indminTOTE = find(TOTE==min(TOTE));
    optLAM = lamarr(indminTOTE,:);  
    totalerr = minTOTE
    xfn = strcat(fnpath, 'CVE_', abname, '_', regtype, '.csv');
    csvwrite(xfn, CVE);
    
   
%     figcv = graphlasso(A, B, XV(indminTOTE,:),INDSV(indminTOTE,:), ...
%                         XT(indminTOTE,:),INDST(indminTOTE,:),optLAM, ...
%                         NONZERO(indminTOTE), abname, 'SAT CV' ,0)
    
lamparam = optLAM;

end

