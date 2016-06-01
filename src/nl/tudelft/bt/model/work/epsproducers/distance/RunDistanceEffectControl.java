package nl.tudelft.bt.model.work.epsproducers.distance;


import nl.tudelft.bt.model.Model;

public class RunDistanceEffectControl extends EpsProducersDistanceEffect {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int seed = 100;
		// parse input
		if (args.length > 0) {
			runWithGraphics = true;
			if (args.length != 2) {
				throw new RuntimeException("Program arguments missing: "
						+ "2 program arguments should be supplied\n"
						+ "1 - output directory\n"
						+ "2 - random number generator seed\n");
			}
			outputDirectory = args[0];
			seed = Integer.parseInt(args[1]);
		} else {
			runWithGraphics = true;
		}
		fEPS_WT = 0.5f;
		specificMassEPS_WT = specificMassX / 6f;
		// these are set for a EPS- mutant
		fEPS_EpsMinus = 1e-6f;
		specificMassEPS_EpsMinus = specificMassEPS_WT;
		//
		timeForWritingToDisk = finitshITerationTime;
		initialParticleNumberWT = 0;
		initialParticleNumberEpsMinus = 61;
		Model.model().setSeed(seed);
		//
		run();
	}
}
