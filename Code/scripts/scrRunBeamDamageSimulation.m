% __scrRunBeamDamageSimulation____ 
% run beam damage simulation and show the loss and gain percentages
% to view simulation results run 
% BeamDamageSimulationViewer(bSimulation.results.resultStruct(n,m))
% with n and m representing the simulation indices

bParams     = BeamDamageParams;
bSimulation = BeamDamageSimulation(bParams);
bSimulation.Run;
% Analyze results 
bdr = BeamDamageResultsAnalysis(bSimulation.results);

