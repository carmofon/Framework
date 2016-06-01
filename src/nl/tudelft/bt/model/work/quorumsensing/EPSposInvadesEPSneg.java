package nl.tudelft.bt.model.work.quorumsensing;

import nl.tudelft.bt.model.Model;

public class EPSposInvadesEPSneg extends EPSposVsEPSneg {

	public static void main(String[] args) {
		if (args.length < 5)
			throw new RuntimeException("First argument: output directory;"
					+ " Second argument: graphics on/off;"
					+ " Third argument: Investment into EPS for EPS+ [0,1];"
					+ " Fourth argument: Investment into EPS for EPS- [0,1];"
					+ " Fifth argument: seed");
		if (!args[0].contains("results"))
			throw new RuntimeException(
					"Output directory doesn't contain the word 'results'");
		outputDirectory = args[0];
		graphics = (Integer.parseInt(args[1]) == 0 ? false : true);
		f_EPSpos = Float.parseFloat(args[2]);
		f_EPSneg = Float.parseFloat(args[3]);
		// seting the random number generator seed
		int seed = Integer.parseInt(args[4]);
		Model.model().setSeed(seed);
		outputEvery = 339;
		run();
	}
}
