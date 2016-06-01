package nl.tudelft.bt.model.work.quorumsensing;

import nl.tudelft.bt.model.Model;



/**
 *  Simulations with EPS neg only to quantify invasiveness
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Nov 8, 2007
 */
public class EPSnegForComparison extends EPSposVsEPSneg {

	public static void main(String[] args) {
		if (args.length < 4)
			throw new RuntimeException("First argument: output directory;"
					+ " Second argument: graphics on/off;"
					+ " third argument: Investment into EPS for EPS- [0,1];"
					+ " fourth argument: seed");
		if (!args[0].contains("results"))
			throw new RuntimeException(
					"Output directory doesn't contain the word 'results'");
		outputDirectory = args[0];
		graphics = (Integer.parseInt(args[1]) == 0 ? false : true);
		f_EPSneg = Float.parseFloat(args[2]);
		f_EPSpos = f_EPSneg; // this strain is not used in computation
		initialParticleNumber_EPSpos = 0;
		initialParticleNumber_EPSneg = 100;
		// seting the random number generator seed
		int seed = Integer.parseInt(args[3]);
		Model.model().setSeed(seed);
		outputEvery = 339;
		run();
	}
}
