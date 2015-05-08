%Multiplierless rotator
%implementation of low complexity multiplierless rotators for MCR
%
%Authors- Shashwat Sanghavi,Rahul Patel,Varsha Teckchandani,
%         Anmol Anubhai,Nishanshi Shukla,Rajan Mehta,Keyur Patel,Raj Jhala
%Group 1
clc;
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=6; %Maximum number of adders allowed
Emax=0.05; %max error
alpha(1)=14; %first angle
alpha(2)=38; %second angle
b=6; %coofficent's word bits

%creating the lookup table which will be used for calculating radius with
%2^i
%In the case of multiple angle rotator, we take those coefficients whose
%radius difference is less that 2*Emax.
%We need to scale this error according to the radius of coeefficient,
%so for finding 2's nearest power we are using lookup table as we are
%using unity scalling.

for i=1:b
    lookup(i)=2^(i);
end

% calling SCR function for both angles,

% Each optimization problem of MCR can be treated as a collection of SCR
% problems

% CORDIC_project refers to a SCR function which takes input of maximum nuber
% of adders, word length, maximum error and an angle and gives output as 
% [X coefficient, Y coefficient, list of all possible coeefficients]

[X,Y,R1]=SCR(N,b,Emax,alpha(1)); 
[X,Y,R2]=SCR(N,b,Emax,alpha(2));

%identifying number of rows in R1 and R2
[h,l]=size(R1); 
[h1,l1]=size(R2);

k=1; % This counter is used for maintaining the pairs of cofficients which lie 
     % within the radius boundry (i.e error criteria explained on line no. 21-26)

% identifying pairs withing similar range beside 2^i arc
for i=1:h
    [A, I] = min(lookup - abs(R1(i))); %identifying the index from the lookup table which contains 2's closest power
    
    for j=1:h1
        Rm=2*Emax*lookup(I);   %scalling the error
        diff = abs(R1(i,2)-R2(j,2)); %calculating the difference between two coefficient's radii
        error = max(R1(i,3),R2(j,3)); %max of error of two coefficients will be the error for that pair
        if Rm > diff  %if the difference between two coefficients is less than scaled error than cosider that pair of coefficients valid
            final(k,:)=[R1(i,1) R2(j,1) error];
            k=k+1;
        end        
    end
end
%identifying the coefficient pair with minimum error from the list of valid
%pairs
[A, I] = min(final(:,3));
%storing the identified pair
coeff(1,:) = [real(final(I,1)) imag(final(I,1)) real(final(I,2)) imag(final(I,2))];

%taking pair which has atleast one odd coefficient out of four
flag = 1;
while(1)
    if(mod(coeff(1,1),2)==0 & mod(coeff(1,2),2)==0 & mod(coeff(1,3),2)==0 & mod(coeff(1,4),2)==0)
        coeff(1,1) = coeff(1,1)/2;
        coeff(1,2) = coeff(1,2)/2;
        coeff(1,3) = coeff(1,3)/2;
        coeff(1,4) = coeff(1,4)/2;
    else
        break; 
    end    
end

coeff %X1,Y1,X2,Y2