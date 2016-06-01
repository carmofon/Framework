package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * @author xavier
 */
public class Kd1minus2_3D extends DetachmentGrowth {
	public static void main(String[] args) {
		System.out.println("running Kd1minus2_3D...");
		// create the application
		Kd1minus2_3D app = new Kd1minus2_3D();
		// set dimensionality to 3D
		DetachmentGrowth.geometry = 3;
		// output
		DetachmentGrowth.outputDirectory = "./Kd1minus2_3D/";
		// Computation parameters
		DetachmentGrowth.systemSize = 1000; // [um]
		// division radius 6 um
		DetachmentGrowth.relativeMaximumRadius = 0.006f;
		//minimum radius
		DetachmentGrowth.relativeMinimumRadius = 
			DetachmentGrowth.relativeMaximumRadius * 0.0001f;
		// boundary layer 200 um
		DetachmentGrowth.relativeBoundaryLayer = 0.2f;
		// other model parameters
		DetachmentGrowth.gridSide = 33; // multigrid grid side
		// kshov is kept
		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]
		// adjust the density
		DetachmentGrowth.densityH = 273.3f;
		// detachment rate coefficient
		DetachmentGrowth.rdetach = 1e-2f; // 300 um 3e-3
		// inoculation
		DetachmentGrowth.initialCellNumber = 2500;
		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(4380); // 1/2 year
		Model.model().setCompulsoryTimeStep(20.0f); // every 20 h output
		// visuals off
		app.run(false);
		System.out.println("Finished Kd1minus2_3D!");
	}
}