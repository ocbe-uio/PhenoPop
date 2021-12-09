% Hill function
function [hill] = hillfunc(b,x) %y = b + (1 - b) ./ (1 + (x./E).^n)

hill = b(2)+ (1-b(2))./(1+(x./b(3)).^b(4));
end
