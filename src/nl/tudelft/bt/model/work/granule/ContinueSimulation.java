/*
 * File created originally on Jun 27, 2005
 */
package nl.tudelft.bt.model.work.granule;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ModelHandler;

/**
 * @author jxavier
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Generation - Code and Comments
 */
public class ContinueSimulation {
	public static void main(String[] args) {
		if (args.length != 2) {
			throw new RuntimeException("Program arguments missing: "
					+ "2 program arguments should be supplied"
					+ " (1 - the output directory,"
					+ " 2 - the filename of model handler dump file)");
		}
		// parse the input arguments
		String outputDir = args[0];
		String fileToRead = args[1];
		try {
			ModelHandler app = ModelHandler.loadFromDumpFile(fileToRead);
			// set the output directory without cleaning contents
			app.restartStateWriters(outputDir);
			// start running
			try {
				// solve the solute concentrations again
				Model.model().solveWithMultigrid();
				// and restart the iterating
				app.startIterating();
				// compute the solute concentration
				//Model.model().solveWithMultigrid();
				// write to disk
				//app.forceWriteState();
			} catch (Exception e) {
				e.printStackTrace();
				app.forceWriteState();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}
}