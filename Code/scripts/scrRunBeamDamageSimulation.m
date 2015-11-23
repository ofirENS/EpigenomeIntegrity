% scrRunBeamDamageSimulation 
% run beam damage simulation and show the loss and gain percentages

bParams     = BeamDamageParams;
bSimulation = BeamDamageSimulation(bParams);
bSimulation.Run;
% Analyze results 
bdr = BeamDamageResultsAnalysis(bSimulation.results);

