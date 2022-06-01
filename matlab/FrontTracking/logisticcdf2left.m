function y = logisticcdf2left(x,a,b,mu,I)
y=a+b.*(1./(1+exp(+(x-mu)/I)).^2);
end
