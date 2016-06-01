package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * 5 replicate runs of the kdet = 1e-2 case
 * 
 * @author xavier
 */
public class Comp1LowDetKd1minus3Replicates extends DetachmentGrowth {
	private static void initiateSystem(int n) {
		int seed = (int) System.currentTimeMillis();
		System.out.println("replicate " + n + "; seed used: " + seed);
		Model.model().setSeed(seed);
		// setting the seed value must come before than reseting the model
		Model.model().reset();
		DetachmentGrowth.outputDirectory = "/Users/jxavier/results/detachmentVsGrowth/"
				+ "comparison1/"
				+ "Comparison1_highdetachment_kd1e-3qmax9.52e-1_replicate"
				+ n
				+ "/";
		// Computation parameters
		DetachmentGrowth.systemSize = 4000; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.0015f;

		DetachmentGrowth.relativeMinimumRadius = DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.05f;

		// other model parameters
		DetachmentGrowth.gridSide = 129; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 1e-3f;// 300 um 3e-3

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 1000;

		//set iteration finish time in h and the iteration time step
		//Model.model().setFinishIterationTime(4000);
		Model.model().setFinishIterationTime(10);
		Model.model().setCompulsoryTimeStep(50.0f);
	}

	public static void main(String[] args) {
		for (int i = 1; i < 6; i++) {
			initiateSystem(i);
			Comp1LowDetKd1minus3Replicates app = new Comp1LowDetKd1minus3Replicates();
			app.run(true);
		}
	}
}