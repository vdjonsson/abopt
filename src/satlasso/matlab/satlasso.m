 function [y] = satlasso(Agood, Bgood, Asat, Bsat, lamarr, regtype)
    % This is where teh debug went 
    % This function solves the augmented lasso problem and returns
    % y =  the refined resgressors
    % x =  non refined regressors,
    % given inputs: 
    % Agood = the sequence corressponding to the "good" data 
    % Bgood = the neutralization values corresponding to the "good" data
    % Asat  = the sequence corresponding to the saturated data 
    % Bsat  = the neutralization values corresponding to the saturated data 
    % lamarr = the array of lambdas to use 

    AT = Agood;
    BT = Bgood;
    
    ATS = Asat;
    BTS = Bsat;
    
    % concatenate good and saturated data to graph later
    ATT = [AT; ATS];

    n = size(AT,2); % col space, max size of regressors

    lam1 = lamarr(1);
    lam2 = lamarr(2);
    lam3 = lamarr(3);
   
    if strcmp(regtype,'elasticnet') == true
         lam4 = lamarr(4);  
    end

    msize = 1/length(BT); % size of the good data   
    
    if strcmp(regtype, 'lasso-5-fold') == true 
        cvx_clear
        cvx_quiet(true) 
        cvx_begin
            variable x(n)
            f = msize*lam1*norm(AT*x-BT, 2) ...
            + msize*lam2*norm(x,1) ...
            + lam3*max([-(ATS*x-BTS);0]);
            minimize(f);
        cvx_end

        for i=1:length(x)
            if abs(x(i)) <= 1e-4
                x(i) = 0;
            end
        end 

        % Refinement step 
        cvx_quiet(true) 
        cvx_begin
            variable y(n)
            f =   lam1*msize*norm(AT*y-BT, 2) ...
            + lam3*max([-(ATS*y-BTS);0]);  
            minimize(f);
            subject to 
            y.*[x==0]==0;
        cvx_end
        
    elseif strcmp(regtype, 'squarelasso') == true 
        cvx_clear
        cvx_quiet(true) 
        cvx_begin
            variable x(n)
            f = msize*lam1*(AT*x-BT)'*(AT*x-BT) ...
                + msize*lam2*norm(x,1) ...
                + lam3*max([-(ATS*x-BTS);0]);
            minimize(f);
        cvx_end

        for i=1:length(x)
            if abs(x(i)) <= 1e-4
                x(i) = 0;
            end
        end 

        % Refinement step 
        cvx_quiet(true) 
        cvx_begin
            variable y(n)
            f =   lam1*msize*(AT*x-BT)'*(AT*x-BT) ...
            + lam3*max([-(ATS*y-BTS);0]);  
            minimize(f);
            subject to 
            y.*[x==0]==0;
        cvx_end
    elseif strcmp(regtype, 'elasticnet') == true 
        cvx_clear
        cvx_quiet(true) 
        cvx_begin
            variable x(n)
            f = msize*lam1*norm(AT*x-BT, 2) ...
            + msize+lam4*x'*x ...
            + msize*lam2*norm(x,1) ...
            + lam3*max([-(ATS*x-BTS);0]);
            minimize(f);
        cvx_end

        for i=1:length(x)
            if abs(x(i)) <= 1e-4
                x(i) = 0;
            end
        end 
        
        y = x;        
        
    end

    for i=1:length(y)
        if abs(y(i)) <= 1e-4
            y(i) = 0;
        end
    end 
end

