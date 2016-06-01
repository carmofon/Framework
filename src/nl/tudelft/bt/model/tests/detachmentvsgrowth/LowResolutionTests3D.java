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
public class LowResolutionTests3D extends DetachmentGrowth {
	public static void main(String[] args) {
		DetachmentGrowth.geometry = 3;
		DetachmentGrowth.densityH = 250;

		DetachmentGrowth.outputDirectory =
			"E:/results/detachmentVsGrowth/" + "lixo/";
//		"E:/results/detachmentVsGrowth/" + "3d/comp1_kd3e-2/";
		// Computation parameters
		DetachmentGrowth.systemSize = 500; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.012f;

		DetachmentGrowth.relativeMinimumRadius =
			DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.1f;

		// other model parameters
		DetachmentGrowth.gridSide = 17; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 1.5e-4f;

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 2500; //2500

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(4000);
		Model.model().setCompulsoryTimeStep(20f);
		//
		DetachmentGrowth.run(args);
	}
}