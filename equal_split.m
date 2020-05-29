function a =  equal_split(x,n) 
    %% a = equal_split(x,n)
    %
    % This function gives the most equal split of the
    % number x into n parts e.g. equal_split(10,3) = [4,3,3]
    %
    % x: integer to be split
    % n: number of parts x is split into.
    % a: vector of length=n, values are in descending order, as equal as possible
    % and sum to x.
    %%
    %If we cannot split the  
    %number into exactly 'N' parts 
    if n == 0
        warning('n is 0, setting to 1');
        n=1;
    end
    
    a = zeros(1,n);
    if(x<n) 
       error('x is less than n');
       
  
    %If x % n == 0 then the minimum  
    %difference is 0 and all  
    %numbers are x / n 
    elseif (mod(x,n) == 0)
        for i=1:n
            a(i)=idivide(x,n);
        end
    else 
        %upto n-(x % n) the values  
        %will be x / n  
        %after that the values  
        %will be x / n + 1 
        zp = n - mod(x,n); 
        pp = idivide(x,n);
        for i=1:n
            if(i> zp)
                a(i)=pp+1;
            else 
                a(i)=pp;
            end
        end
    end    
    a = sort(a,'descend');
end