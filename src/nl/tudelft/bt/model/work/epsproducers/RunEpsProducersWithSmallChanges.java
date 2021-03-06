package nl.tudelft.bt.model.work.epsproducers;

import nl.tudelft.bt.model.exceptions.ModelRuntimeException;

/**
 * Runs EpsProducers simulation where WT competes with EPS- mutant with variable
 * initial frequency of WT
 * 
 * @author jxavier
 * 
 */
public class RunEpsProducersWithSmallChanges {
	private static final float SMALLCHANGE = 0.1f;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		if (args.length != 8) {
			throw new RuntimeException("Program arguments missing: "
					+ "8 program arguments should be supplied\n"
					+ "1 - output directory"
					+ "2 - flag for running with graphics\n"
					+ "3 - fEPS of WT\n"
					+ "4 - desity of biomass / density of EPS in WT\n"
					+ "5 - original frequency of WT\n"
					+ "6 - oxygen bulk concentration\n"
					+ "7 - boundary layer thickness\n"
					+ "8 - number of replicate tag\n" + ")");
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
		EpsProducers.fEPS_EpsMinus = EpsProducers.fEPS_WT - SMALLCHANGE;
		// determine the original frequencies
		// parse the initial frequencies
		int total = 100;
		EpsProducers.initialParticleNumberWT = Math.round(total
				* Float.parseFloat(args[4]));
		EpsProducers.initialParticleNumberEpsMinus = total
				- EpsProducers.initialParticleNumberWT;
		//
		EpsProducers.specificMassEPS_EpsMinus = EpsProducers.specificMassEPS_WT;
		EpsProducers.oxygenBulkConcentration = Float.parseFloat(args[5]) * 1e-3f;
		EpsProducers.relativeBoundaryLayer = Float.parseFloat(args[6])
				/ EpsProducers.systemSize;
		// check if replicate flag is in correct format
		try {
			Integer.parseInt(args[7]);
		} catch (NumberFormatException e2) {
			throw new ModelRuntimeException("Replicate number " + args[7]
					+ " not in correct format");
		}
		EpsProducers.run();
	}
}
