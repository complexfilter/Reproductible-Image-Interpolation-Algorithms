% 2X enlarge

function y=sri(x,level)

   for l=1:level
   y1=sri1(x);y=sri2(y1);
   x=y;
   end

