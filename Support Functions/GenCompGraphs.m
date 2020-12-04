function [Block, Comp] = GenCompGraphs(GRAPH)

Block = {};
Comp  = [];

for i = 1:numel(GRAPH)
    switch GRAPH(i).Type
        case 'Tank'
            Model = NonpressurizedLiquidTank;
            Model.Name = GRAPH(i).Name;
            Model.fluid_props = 'JP8';
            Model.initial_mass_lb = 100;
            Model.initial_volume_in3 = 0;
            Model.initial_temperature_R = 580;
            Model.tank_volume_in3 = 600000;
            skip = 0;
        case 'Heat Exchanger'
            
            % Build Heat Exchanger
            Model = HeatExchanger;
            Model.Name = GRAPH(i).Name;
            Model.fluid_props = 'JP8';
            Model.HTC = 100;
            skip = 0;
            
        case 'Heat Source'
            
            % Build External Heat Graph
            Model = ExternalHeat;
            Model.Name = GRAPH(i).Name;
            Model.fluid_props = 'JP8';
            skip = 0;
            
        case 'Junction'
            
            % Build Junction
            Model = MergeManifold;
            Model.Name = GRAPH(i).Name;
            Model.fluid_props = 'JP8';
            Model.N_in = length(GRAPH(i).UpVertex);
            Model.N_out = length(GRAPH(i).DownVertex);
            skip = 0;
            
        case 'Split'
            
            % Build Junction
            Model = MergeManifold;
            Model.Name = GRAPH(i).Name;
            Model.fluid_props = 'JP8';
            Model.N_in = length(GRAPH(i).UpVertex);
            Model.N_out = length(GRAPH(i).DownVertex);
            skip = 0;
            
        otherwise
            skip = 1;
            
    end
    
    if skip == 0
        Model.generateGraph; % call this function to generate a component graph
        
        Block = [Block {Model}];
        Comp  = [Comp Model.Comp_Graph];
    end
end

figure
for i = 1:numel(Comp)
    subplot(ceil(numel(Comp)/2),2,i)
    G = digraph(Comp(i).E(:,1),Comp(i).E(:,2));
    h = plot(G);
    labeledge(h,Comp(i).E(:,1)',Comp(i).E(:,2)',1:Comp(i).Ne);
    title(Block{i}.Name)
end
    
    
end