package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * 5 replicate runs of the kdet = 1e-2 case
 * 
 * @author xavier
 */
public class Comp2VariableKdet3DRuns extends DetachmentGrowth {
	private static void initiateSystem(float kdet) {
		// reset the model
		Model.model().reset();
		//
		int seed = (int) System.currentTimeMillis();
		System.out.println("kdet " + kdet + ";");
//		DetachmentGrowth.outputDirectory = "D:/JOAO/results/"
//			+ "comparison2/"
//			+ "kd"+ kdet +"qmax9.52e-1/";
		DetachmentGrowth.outputDirectory = "E:/results/detachmentVsGrowth/3d/"
			+ "kd"+ kdet +"qmax9.52e-1/";
		// Computation parameters
		DetachmentGrowth.geometry = 3;
		
		DetachmentGrowth.systemSize = 500; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.012f;

		DetachmentGrowth.relativeMinimumRadius =
			DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.4f;

		// other model parameters
		DetachmentGrowth.gridSide = 17; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = kdet;// 300 um 3e-3

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 2500; //2500

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(4000);
		Model.model().setCompulsoryTimeStep(20f);

	}

	public static void main(String[] args) {
		//float[] kdets = {1.5e-6f, 1.5e-5f, 1.5e-4f};
		float[] kdets = {1.5e-4f};
		for (int i = 0; i < kdets.length; i++) {
			initiateSystem(kdets[i]);
			DetachmentGrowth.run(args);
		}
	}
}