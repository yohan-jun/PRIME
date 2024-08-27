function [g] = slewRate_smooth(g_input)


g = g_input ; % Gauss/cm

ng = abs(g./max(g));
indx = find((ng) >=0.99);
g(1:indx(1)) = linspace(g(1),g(indx(1)),length(1:indx(1)));
iL = length(indx);
g(indx(iL):length(g))=linspace(g(iL),g(end),length(indx(iL):length(g)));