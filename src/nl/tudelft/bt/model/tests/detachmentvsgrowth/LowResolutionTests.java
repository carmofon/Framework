/*
 * Created on Jun 30, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * @author xavier
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 */
public class LowResolutionTests extends DetachmentGrowth {
	public static void main(String[] args) {
		DetachmentGrowth.outputDirectory =
			"/Users/jxavier/results/detachmentVsGrowth";
		
		// uptake rate
		DetachmentGrowth.qMax = 0.8f; //[gCOD_S/gCOD_X/h] Beun 
		// Computation parameters
		DetachmentGrowth.systemSize = 1000; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.006f;

		DetachmentGrowth.relativeMinimumRadius =
			DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.2f;

		// other model parameters
		DetachmentGrowth.gridSide = 33; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 1.5e-5f;

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 250;

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(8760);
		Model.model().setCompulsoryTimeStep(50.0f);
		//
		LowResolutionTests app = new LowResolutionTests();
		app.run(true);
	}
}