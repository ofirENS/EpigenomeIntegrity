% scrRunBeamDamageSimulation 
bParams     = BeamDamageParams;
bSimulation = BeamDamageSimulation(bParams);
bSimulation.Run;
% Analyze results 
bdr = BeamDamageResultsAnalysis(bSimulation.results);

