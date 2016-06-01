package nl.tudelft.bt.model.work.epsproducers;

import nl.tudelft.bt.model.exceptions.ModelRuntimeException;

/**
 * Runs EpsProducers simulation where WT competes with EPS- mutant with variable
 * initial frequency of WT
 * 
 * @author jxavier
 * 
 */
public class RunEpsProducersWithFrequencies {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		EpsProducers.runWithGraphics = false;
		if (args.length != 6) {
			throw new RuntimeException("Program arguments missing: "
					+ "4 program arguments should be supplied\n"
					+ "1 - output directory"
					+ "2 - flag for running with graphics\n"
					+ "3 - fEPS of WT\n"
					+ "4 - desity of biomass / density of EPS in WT\n"
					+ "5 - initial frequency of WT in [0,1] range\n"
					+ "6 - number of replicate tag\n" + ")");
		}
		// pass output directory
		EpsProducers.outputDirectory = args[0];
		// parse GUI flag
		int arg1 = Integer.parseInt(args[1]);
		switch (arg1) {
		case 0:
			EpsProducers.runWithGraphics = false;
			break;
		case 1:
			EpsProducers.runWithGraphics = true;
			EpsProducers.timeForWritingToDisk = 1.0f;
			break;
		default:
			throw new RuntimeException("second program" + " argument must be 0"
					+ " (for running with no graphics) "
					+ "or 1 (for running with graphics)");
		}
		// parse model parameters
		EpsProducers.fEPS_WT = Float.parseFloat(args[2]);
		EpsProducers.specificMassEPS_WT = EpsProducers.specificMassX
				/ Float.parseFloat(args[3]);
		// these are set for a EPS- mutant
		EpsProducers.fEPS_EpsMinus = 1e-6f;
		EpsProducers.specificMassEPS_EpsMinus = EpsProducers.specificMassEPS_WT;
		// parse the initial frequencies
		int total = 100;
		EpsProducers.initialParticleNumberWT = Math.round(total
				* Float.parseFloat(args[4]));
		EpsProducers.initialParticleNumberEpsMinus = total
				- EpsProducers.initialParticleNumberWT;
		// check if replicate flag is in correct format
		try {
			Integer.parseInt(args[5]);
		} catch (NumberFormatException e2) {
			throw new ModelRuntimeException("Replicate number " + args[8]
					+ " not in correct format");
		}
		EpsProducers.run();
	}
}
