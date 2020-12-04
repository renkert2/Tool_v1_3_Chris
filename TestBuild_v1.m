clear all
close all
clc

set(0,'defaultAxesFontSize' , 12)
set(0,'defaultTextFontSize' , 12)

%% Information Provided from simulink Model
load ToySystem
plot(G)

%% Build Component Graphs
[Block, Comp] = GenCompGraphs(cGRAPH); 

%% Get Graph Interconnections
ConnectE = ExtractExConn(Comp,cGRAPH); 
ConnectV = ExtractVxConn(Comp,ConnectE); 

%% Build the System Graph
Sys = GenSysGraph(Comp,ConnectV,ConnectE);

%%
Sys = SymbolicSolver(Sys);

% test

%% Linearize the System Graph
% x0 = [60; 70; 580; 580; 585; 590; 595; 600; 605; 0; 0; 610; 615]; 
% u0 = [.1 .2 .4 .1 .2 .3]';
% Sys = LinearizeGraph(Sys,x0,u0);










%% To Do

