package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * 5 replicate runs of the kdet = 1e-2 case
 * 
 * @author xavier
 */
public class Comp2ChangingCO extends DetachmentGrowth {
	private static void initiateSystem(float co) {
		System.out.println("New run: CO " + co);
		// setting the seed value must come before than reseting the model
		Model.model().reset();

		// uptake rate
		DetachmentGrowth.oxygenBulkConcentration = co; //[gO]
		DetachmentGrowth.outputDirectory = "D:/JOAO/results/"
				+ "Comparison2_kd3e-4qmax0.952Co"+ co +"/";

		// Computation parameters
		DetachmentGrowth.systemSize = 4000; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.0015f;

		DetachmentGrowth.relativeMinimumRadius = DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.05f;

		// other model parameters
		DetachmentGrowth.gridSide = 129; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 3e-4f;// 300 um 3e-3

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 1000;

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(4000/co);
		//Model.model().setFinishIterationTime(10/qmax);
		Model.model().setCompulsoryTimeStep(50.0f/co);
	}

	public static void main(String[] args) {
		float [] co = {9e-3f};
		for (int i = 0; i < co.length; i++) {
			initiateSystem(co[i]);
			DetachmentGrowth.run(args);
		}
	}
}