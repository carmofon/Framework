package nl.tudelft.bt.model.work.epsproducers;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;

/**
 * Run a simulations without gradients to provide comparison for the cases with
 * gradients
 * 
 * @author jxavier
 * 
 */
public class RunEpsProducersNoGradients {
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// pass output directory
		EpsProducers.outputDirectory = "C:/results/no_gradients/";
		EpsProducers.runWithGraphics = true;
		// parse model parameters
		EpsProducers.fEPS_WT = 0.55f;
		EpsProducers.specificMassEPS_WT = 6f;
		// these are set for a EPS- mutant
		EpsProducers.fEPS_EpsMinus = 1e-6f;
		EpsProducers.specificMassEPS_EpsMinus = EpsProducers.specificMassEPS_WT;
		// parse the initial frequencies
		EpsProducers.initialParticleNumberWT = 50;
		EpsProducers.initialParticleNumberEpsMinus = 50;
		//set iteration time to 1h 
		EpsProducers.timeForWritingToDisk = 0.1f;
		EpsProducers.finitshITerationTime = 8f;
		// set oxygen diffusivity to a value that's icredibly high
		EpsProducers.oxygenDiffusivity *= 1000f;
		// that's it, now run
		EpsProducers.run();
	}
}
