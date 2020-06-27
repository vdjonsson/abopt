function [P] = satlassostabselection(fullA, fullB, satindex, abname, numperturbations, lamarr)

% This function does stability selection on the full data set, including
% saturation points and uses the augmented lasso function that penalizes 
% the saturated data 
% Given the number of data perturbations, partition into training set 
% and validation set.  Each set is picked randomly without replacement, 
% and a total of n/2 samples are picked at each run. 


    A = fullA; 
    B = fullB; 
    bstr = abname; 
    lenlam = size(lamarr,1);
    
    % First separate all the data into good data and saturated data     
    sizealldata = size(fullB,1);
    n = size(fullA,2); %col space, max size of regressors
    
    % This is the probability of occurrence of each regressor per lambda
    % value 
    P = [];
    
    for jj=1:lenlam
        
        Plam = [];
        
        lamarr(jj,:)
        
        jj;
        lam1 = lamarr(jj, 1);
        lam2 = lamarr(jj, 2);
        lam3 = lamarr(jj, 3);
    
        % Indices for good and saturated data 
        BTi = 1:satindex;
        BTSi = satindex+1:sizealldata;
        
        satval = fullB(satindex+1);
       
        pctrain = 1/2;
     
        sizesatdata = length(BTSi);
        sizedata = sizealldata-sizesatdata; 
        
        trainsize = floor(sizealldata*pctrain); 
        
        
        % For each training partition pick off x% points from good data 
        % and y% points from the saturated data, where x and y are
        % calculated depending on the ratio of good and bad data 
        
        pcsatincl = sizesatdata/sizealldata; 
         
        trainsatsize = floor(trainsize*pcsatincl); 
        traingoodsize = trainsize - trainsatsize;
        
        if traingoodsize > length(BTi)
            traingoodsize = traingoodsize -1;
        end
     
        Pk = zeros(n,1);
        for j=1:numperturbations
            
            trgoodinds = datasample(BTi,traingoodsize,'Replace', false); 
            trsatinds = datasample(BTSi,trainsatsize);  
            
            % Now get validation indices different from training ones
            i= 1; 
            
            trallinds = [trgoodinds, trsatinds];
            
       
            
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
    
            [y] = satlasso(AT, BT, ATS, BTS, lamarr(jj,:));
            
            % Calculate the occurence of each regressor 
            Pik = abs(y)>0;
            Pk = Pk + Pik;      
        end   
        Pk = Pk./numperturbations;
   
    P = [P Pk];
    
    end   
end

