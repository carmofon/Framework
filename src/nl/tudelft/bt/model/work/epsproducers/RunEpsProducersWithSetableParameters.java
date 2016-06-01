package nl.tudelft.bt.model.work.epsproducers;

import nl.tudelft.bt.model.exceptions.ModelRuntimeException;

/**
 * Runs EpsProducers simulation where WT competes with EPS- mutant with variable
 * initial frequency of WT
 * 
 * @author jxavier
 * 
 */
public class RunEpsProducersWithSetableParameters {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		if (args.length != 9) {
			throw new RuntimeException("Program arguments missing: "
					+ "4 program arguments should be supplied\n"
					+ "1 - output directory"
					+ "2 - flag for running with graphics\n"
					+ "3 - fEPS of WT\n"
					+ "4 - desity of biomass / density of EPS in WT\n"
					+ "5 - fEPS in EPS- mutant\n"
					+ "6 - desity of biomass / density of EPS in EPS- mutant\n"
					+ "7 - oxygen bulk concentration\n"
					+ "8 - boundary layer thickness\n"
					+ "9 - number of replicate tag\n" + ")");
		}
		// pass output directory
		EpsProducers.outputDirectory = args[0];
		// check if replicate flag is in correct format
		try {
			Integer.parseInt(args[8]);
		} catch (NumberFormatException e2) {
			throw new ModelRuntimeException("Replicate number " + args[8]
					+ " not in correct format");
		}
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
		EpsProducers.fEPS_EpsMinus = Float.parseFloat(args[4]);
		EpsProducers.specificMassEPS_EpsMinus = EpsProducers.specificMassX
				/ Float.parseFloat(args[5]);
		EpsProducers.oxygenBulkConcentration = Float.parseFloat(args[6]) * 1e-3f;
		EpsProducers.relativeBoundaryLayer = Float.parseFloat(args[7])
				/ EpsProducers.systemSize;
		//
		EpsProducers.run();
	}
}
