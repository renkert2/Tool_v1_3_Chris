function [x] = simA(in,Sys)

xd = in(1:4);
xe = in(5:8);
u  = in(9:14);
ee = in(15:17);

x = Sys.solveAlg(xd,xe,u,ee);


end

