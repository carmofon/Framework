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
 */
public class Comp1HighDetKd3eminus4 extends DetachmentGrowth {
	public static void main(String[] args) {
		DetachmentGrowth.outputDirectory = "E:/results/detachmentVsGrowth/" +
				"comparison1/Comparison1_lowdetachment_kd3e-4qmax9.52e-1/";
		// Computation parameters
		DetachmentGrowth.systemSize = 4000; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.0015f;

		DetachmentGrowth.relativeMinimumRadius = 
			DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.05f;

		// other model parameters
		DetachmentGrowth.gridSide = 129; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 3e-4f;// 300 um 3e-3

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 1000;

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(8760);
		Model.model().setCompulsoryTimeStep(50.0f);
		//
		DetachmentGrowth.run(args);
	}
}