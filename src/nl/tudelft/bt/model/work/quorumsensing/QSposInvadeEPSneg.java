package nl.tudelft.bt.model.work.quorumsensing;

import nl.tudelft.bt.model.Model;

public class QSposInvadeEPSneg extends QSposVsEPSpos {

	public static void main (String[] args) {
		if (args.length < 5) throw new RuntimeException("First argument: output directory; Second argument: graphics on/off " +
				"Third argument: AI production rate [0,1]; Fourth argument: Initial # of QSpos (out of 100)" +
				"Fifth argument: seed");
		if (!args[0].contains("results")) throw new RuntimeException("Output directory doesn't contain the word 'results'");
		outputDirectory = args[0];
		graphics = (Integer.parseInt(args[1]) == 0 ? false : true);
		float alpha_tilde = Float.parseFloat(args[2]);
		float L2 = 100*100;
		QSthreshold = alpha_tilde*L2*aiProductionRate*specificMassX/autoinducerDiffusivity;
		initialParticleNumber_QSpos = Integer.parseInt(args[3]);
		initialParticleNumber_QSneg = 100 - initialParticleNumber_QSpos;
		f_QSneg = 0f;
		initialBiomassQSpos = 1;
		initialBiomassQSneg = 1;
		// seting the random number generator seed
		int seed = Integer.parseInt(args[4]);
		Model.model().setSeed(seed);
		outputEvery = 300;
		run();
	}
}
