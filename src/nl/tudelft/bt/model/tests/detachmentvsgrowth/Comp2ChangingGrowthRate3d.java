package nl.tudelft.bt.model.tests.detachmentvsgrowth;

import nl.tudelft.bt.model.Model;

/**
 * 5 replicate runs of the kdet = 1e-2 case
 * 
 * @author xavier
 */
public class Comp2ChangingGrowthRate3d extends DetachmentGrowth {
	private static void initiateSystem(float qmax) {
		System.out.println("New run: max " + qmax);
		// setting the seed value must come before than reseting the model
		Model.model().reset();

		// uptake rate
		DetachmentGrowth.qMax = qmax; //[gCOD_S/gCOD_X/h] Beun
		DetachmentGrowth.outputDirectory = "E:/results/detachmentVsGrowth/"
				+ "comparison2VariableGrowthRates/" + "kd3e-3qmax"
				+ DetachmentGrowth.qMax + "_3D/";

		// Computation parameters
		DetachmentGrowth.geometry = 3;
		
		DetachmentGrowth.systemSize = 600; // [um]

		DetachmentGrowth.relativeMaximumRadius = 0.0114f;

		DetachmentGrowth.relativeMinimumRadius = DetachmentGrowth.relativeMaximumRadius * 0.0001f;

		DetachmentGrowth.relativeBoundaryLayer = 0.333f;

		// other model parameters
		DetachmentGrowth.gridSide = 17; // multigrid grid side

		DetachmentGrowth.kShov = 1.0f; // shoving parameter[dim/less]

		DetachmentGrowth.rdetach = 3e-3f;// 300 um 3e-3

		// detachment constant[g/l/h]
		DetachmentGrowth.initialCellNumber = 1500;

		//set iteration finish time in h and the iteration time step
		Model.model().setFinishIterationTime(4000);
		//Model.model().setFinishIterationTime(10 / qmax);
		Model.model().setCompulsoryTimeStep(50.0f / qmax);
	}

	public static void main(String[] args) {
		float[] qmaxes = { 0.1f, 0.5f, 0.952f};
		for (int i = 0; i < qmaxes.length; i++) {
			initiateSystem(qmaxes[i]);
			DetachmentGrowth.run(args);
		}
	}
}