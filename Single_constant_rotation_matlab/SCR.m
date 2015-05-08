%Multiplierless rotator
%implementation of low complexity multiplierless rotators for SCR
%
%Authors- Shashwat Sanghavi,Rahul Patel,Varsha Teckchandani,
%         Anmol Anubhai,Nishanshi Shukla,Rajan Mehta,Keyur Patel,Raj Jhala
%Group 1

%Defining Inputs
%where numberOfAdders is maximum number of adders allowed, b is word
%length, Emax is the maximum error allowed, alpha1 is the input angle
%provided
function [X,Y,finalR]=SCR(numberOfAdders,b,Emax,alpha1)
    m=tand(alpha1);   %slope of the line at angle alpha1
    Xline=0:2^(b-1); %initializing the design space consisting of all possible finite word length values
    Yline=m*Xline;   %initializing for a kernel of m angles.
    x=0:(2^(b-1)-1);
    y=x;
    figure;
    subplot(1,3,1);
    %showing the line with desired angle
    plot(Xline,Yline);
    hold on;

%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Step-1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for a=1:length(x)     %creating the subspace visual
        for b=1:length(y)
            space(a,b)=x(a)+1i*y(b);
            plot(real(space(a,b)),imag(space(a,b)),'.','MarkerSize',12);
            hold on;
        end
    end
 %%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%step-2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta=asind(Emax); %The maximum delta at which the angle can differ
    idx=1;
    subplot(1,3,2);
    %showing the lines which only differ by delta
    plot(Xline,Yline); 
    hold on;

    for a=1:length(x)
        for b=1:length(y)
            % The criteria for a coefficient to be valid is that 
            %the tan inverse of the Sine by Cosine ratio minus alpha should be less than delta
            if (abs(atand(imag(space(a,b))/real(space(a,b)))-alpha1)<delta) 
                interested(idx)=space(a,b); %storing the valid coefficients
                idx=idx+1;
                 plot(real(space(a,b)),imag(space(a,b)),'.','MarkerSize',12); % creating a visual of the selected coefficients
                 hold on;
            end
        end
    end
  %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%step-3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=1;
    moreInt = [];
    for i=b:-1:0
        R_curr = 2^(i);
        R_Emax = R_curr * Emax; %Maximum permissible value for each value of R_curr
        for j=1:length(interested)
             if (abs(R_curr - abs(interested(j)))  <= R_Emax) % Any coef?cient whose scaling differs more 
                                                              %than R_Emax from R_curr is discarded
                moreInt(k)=interested(j); %Storing the coefficients that satisfy the criteria of being less than R_Emax
                k=k+1;
            end
        end    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fixed R%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     R_curr = 10;
    %     R_Emax = R_curr * Emax;
    %     for j=1:length(interested)
    %        %  if (abs(R_curr - abs(interested(j)))  <= R_Emax)
    %          if (abs(interested(j))  <= R_curr+R_Emax)
    % 
    %             moreInt(k)=interested(j);
    %             k=k+1;
    %     end    
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%step-4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %A lookup table of all C and S values along with the number of adders
    %i.e. C is the real part and S is the imaginary part.
    %for each rotation is prepared.
    %Coef?cients which require more than the allowed number of adders are discarded.
    result = [];
    Ec = real(moreInt)./abs(moreInt) - cosd(alpha1); 
    Es = imag(moreInt)./abs(moreInt) - sind(alpha1);
    Error = (Ec.^(2)+Es.^(2)).^(1/2);

    for i=1:length(moreInt)
        AMr=length(strfind(int2str(de2bi(real(moreInt(i)))), '1'))-1;
        AMi=length(strfind(int2str(de2bi(imag(moreInt(i)))), '1'))-1;
        if(real(moreInt(i))==0)
            AR(i)=2+AMi; %If the rotation coef?cient is real,
            %P=C the rotator is reduced to two real multiplications
        elseif(imag(moreInt(i))==0)
            AR(i)=2+AMr;Likewise, % if the coef?cient is a pure imaginary number, 
            %P=jS the rotation has two real constant multiplications and the number of adders is given by the formula 
        else
           AR(i)=2+AMi+AMr; %The shift-and-add implementation depends on the rotation angle. 
           %In general, a rotation of P=C+jS is calculated and the 
           %total number of adders of the shift-and-add implementation is
           %obtained by this formula.
        end
    end
    
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%step 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    table=reshape(moreInt,[length(moreInt),1]); %the SCR result table is prepared
    result=[table,(abs(table)), Error' , AR'];
    ii=1;
    [Hresult,Lresult]=size(result);
    for i=1:Hresult
        if result(i,4)<=numberOfAdders
        finalR(ii,:)=result(i,:); 
        ii=ii+1;
        end
    end
    
    subplot(1,3,3); 
    plot(Xline,Yline); 
    hold on;
    plot(finalR,'.');
    [c,i]=min(finalR(:,3));%.Typically the one with the smallest rotation error is selected 
                           %as the number of adders are within the speci?cation boundaries
    X=real(finalR(i,1));
    Y=imag(finalR(i,1));
end