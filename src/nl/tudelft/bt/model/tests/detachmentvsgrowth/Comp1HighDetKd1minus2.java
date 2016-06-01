package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * @author xavier
 */
public class Comp1HighDetKd1minus2 extends DetachmentGrowth {
	public static void main(String[] args) {
		DetachmentGrowth.outputDirectory = "d:/results/detachmentVsGrowth/"
				+ "comparison1/"
				+ "Comparison1_highdetachment_kd1e-2qmax9.52e-1/";
		// Computation parameters
		DetachmentGrowth.systemSize = 4000; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.0015f;

		DetachmentGrowth.relativeMinimumRadius = DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.05f;

		// other model parameters
		DetachmentGrowth.gridSide = 129; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 1e-2f;// 300 um 3e-3

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 1000;

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(8760);
		Model.model().setCompulsoryTimeStep(50.0f);
		//
		DetachmentGrowth.run(args);
	}
}