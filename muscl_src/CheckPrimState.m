%Checks that primitive state W is legit
%Ie, pos density and pressure
function W = CheckPrimState(W)
global SMALL_NUM;
global Pmin;
global Dmin;



Ind = ( W(1,:,:) < Dmin ) | ( W(4,:,:) < Pmin );
if sum(Ind(:)) > 0
    %keyboard
    W(1,Ind) = Dmin;
    W(4,Ind) = Pmin;
    W(2,Ind) = 0.0; W(3,Ind) = 0.0;
end
Ind = ~isreal(W);
if sum(Ind(:)) > 0
    W = real(W);
end