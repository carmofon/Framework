package nl.tudelft.bt.model.work.epsproducers.mutation;

import javax.print.attribute.standard.Finishings;

import nl.tudelft.bt.model.Model;

public class RunWithMutation extends EpsProducersMutationImpact {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int seed = 100;
		// parse input
		if (args.length > 0) {
			runWithGraphics = false;
			if (args.length != 6) {
				throw new RuntimeException("Program arguments missing: "
						+ "6 program arguments should be supplied\n"
						+ "1 - output directory\n"
						+ "2 - time for mutation [in h]\n"
						+ "3 - radius of particles to follow\n"
						+ "4 - random number generator seed\n"
						+ "5 - flag for mutation [1 - mutation, 0 - control]\n"
						+ "6 - finish iteration time\n");
			}
			outputDirectory = args[0];
			mutationTime = Float.parseFloat(args[1]);
			neighborhoodRaius = Integer.parseInt(args[2]);
			seed = Integer.parseInt(args[3]);
			// parse the flag
			mutationFlag = (Integer.parseInt(args[4]) == 1);
			//
			finitshITerationTime = Float.parseFloat(args[5]);
		} else {
			runWithGraphics = true;
			initialParticleNumberWT = 40;
			timeForWritingToDisk = 5f;
		}
		fEPS_WT = 0.5f;
		specificMassEPS_WT = specificMassX / 6f;
		// these are set for a EPS- mutant
		fEPS_EpsMinus = 1e-6f;
		specificMassEPS_EpsMinus = specificMassEPS_WT;
		// parse the initial frequencies
		int total = 100;
		initialParticleNumberWT = Math.round(total * 0.50f);
		initialParticleNumberEpsMinus = total - initialParticleNumberWT;
		//
		timeForWritingToDisk = finitshITerationTime;
		Model.model().setSeed(seed);
		//

		//
		run();
	}
}
