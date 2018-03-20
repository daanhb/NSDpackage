%adds all NSD45 relevant paths:
% %addpath NSD;
% addpath NSD/pathChoosing;
% addpath NSD/pathFinding;
% addpath NSD/quadrature;
% addpath NSD/rootFinding;
% addpath NSD/singularPoints;
% addpath NSD/DaveCode


addpath pathChoosing;
addpath pathFinding;
addpath quadrature;
addpath rootFinding;
addpath singularPoints;
addpath DaveCode

%some of Daan's old hoipack paths get used too:
% addpath hoipack;
% addpath hoipack/expansion;
% addpath hoipack/quadrules;
% addpath hoipack/OPQ;

%Moved Daan's hoipack stuff into quadrature folder:
addpath quadrature/expansion;
addpath quadrature/quadrules;
addpath quadrature/OPQ;

%(should probably go through these at some point and delete what's not
%needed)