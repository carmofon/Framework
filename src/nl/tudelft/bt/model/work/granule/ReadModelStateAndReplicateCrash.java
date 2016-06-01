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
public class ReadModelStateAndReplicateCrash {
	private static String fileToRead = "results/UpdatedAquasimModel/modelhandler_iteration00031842.dump";
	private static String outputDirectory = "results/ReplicateCrash";

	public static void main(String[] args) {
		try {
			ModelHandler app = ModelHandler.loadFromDumpFile(fileToRead);
			app.intializeStateWriters(outputDirectory);
			try {
				// start running
				Model.model().solveWithMultigrid();
				//app.startIterating();
				app.forceWriteState();
			} catch (RuntimeException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
				app.forceWriteState();
			}
			
			// compute the solute concentration
			//Model.model().solveWithMultigrid();
			// write to disk
			//app.forceWriteState();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}
}