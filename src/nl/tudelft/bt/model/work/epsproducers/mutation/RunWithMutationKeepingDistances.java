package nl.tudelft.bt.model.work.epsproducers.mutation;

import javax.print.attribute.standard.Finishings;

import nl.tudelft.bt.model.Model;

public class RunWithMutationKeepingDistances extends
		EpsProducersMutationKeepingDistances {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int seed = 100;
		// parse input
		if (args.length > 0) {
			runWithGraphics = false;
			if (args.length != 5) {
				throw new RuntimeException("Program arguments missing: "
						+ "6 program arguments should be supplied\n"
						+ "1 - output directory\n"
						+ "2 - time for mutation [in h]\n"
						+ "3 - random number generator seed\n"
						+ "4 - flag for mutation [1 - mutation, 0 - control]\n"
						+ "5 - finish iteration time\n");
			}
			outputDirectory = args[0];
			mutationTime = Float.parseFloat(args[1]);
			seed = Integer.parseInt(args[2]);
			// parse the flag
			mutationFlag = (Integer.parseInt(args[3]) == 1);
			//
			finitshITerationTime = Float.parseFloat(args[4]);
		} else {
			runWithGraphics = true;
			initialParticleNumberWT = 40;
			timeForWritingToDisk = 5f;
			mutationTime = 0;
		}
		fEPS_WT = 0.5f;
		specificMassEPS_WT = specificMassX / 6f;
		// these are set for a EPS- mutant
		fEPS_EpsMinus = 1e-6f;
		specificMassEPS_EpsMinus = specificMassEPS_WT;
		// parse the initial frequencies
		initialParticleNumberWT = 50;
		initialParticleNumberEpsMinus = 50;
		//
		timeForWritingToDisk = finitshITerationTime;
		Model.model().setSeed(seed);
		//
		//
		run();
	}
}
