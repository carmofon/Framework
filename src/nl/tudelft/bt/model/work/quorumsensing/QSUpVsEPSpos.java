package nl.tudelft.bt.model.work.quorumsensing;

import nl.tudelft.bt.model.Model;

/**
 * Competition betwen a QS strain (upregulation of EPS) and EPS+
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Nov 12, 2007
 */
public class QSUpVsEPSpos extends QSposVsEPSpos {

	public static void main(String[] args) {
		if (args.length < 4)
			throw new RuntimeException("First argument: output directory;"
					+ " Second argument: graphics on/off "
					+ "Third argument: AI production rate [0,1]; Fourth "
					+ " argument: seed");
		if (!args[0].contains("results"))
			throw new RuntimeException(
					"Output directory doesn't contain the word 'results'");
		outputDirectory = args[0];
		graphics = (Integer.parseInt(args[1]) == 0 ? false : true);
		float alpha_tilde = Float.parseFloat(args[2]);
		float L2 = 100 * 100;
		QSthreshold = alpha_tilde * L2 * aiProductionRate * specificMassX
				/ autoinducerDiffusivity;
		// seting the random number generator seed
		int seed = Integer.parseInt(args[3]);
		Model.model().setSeed(seed);
		//IMPORTAT:  sets QS behavior to up-regulation 
		downRegulation = false;
		initialBiomassQSpos = 1;
		initialBiomassQSneg = 1;
		outputEvery = 300;
		run();

	}

}