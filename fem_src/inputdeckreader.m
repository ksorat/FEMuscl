function [NodeList,ElList,DispBC,numnodes,numels] = inputdeckreader(filename)

%% Find Input Deck Sections

[dispbcline,elementline,nodeline,dispbcend,elementend,nodeend] = LineFind(filename); % Find header lines


%% Read Nodes

[NodeList,numnodes] = readnodes(filename,nodeline,nodeend); % Read in nodes

%% Read Displacement BCs

DispBC = readdispbcs(filename,dispbcline,dispbcend); % Read in displacement boundary conditions

%% Read Elements

[ElList,numels] = readelements(filename,elementline,elementend); % Read in elements


