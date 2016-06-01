package nl.tudelft.bt.model.work.tumordemonstration;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.work.relatedness.StepReturnGoodOnly;
import nl.tudelft.bt.model.work.relatedness.StepReturnSubstrateAndGood;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;

/**
 * Simulates growth of two strains with differential public good production in a
 * vertical slice biofilm cross-section
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class ESSGoodsCleanRadial3D extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/careynadell/practice";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 3;

	// These are all now defined in the main function
	protected static boolean gradientsOn;

	protected static float simulationFinishTime; // [h]

	protected static float substrateBulkConcentration;

	// 1. Reaction parameters
	// *********************************************************************************

	// max growth rate on substrate alone
	// private static float uMax = 1.5f;
	private static float uMax = 1f;

	private static float Ks = 3.5e-5f;

	// Effect of public goods upon growth rate
	private static float benefit = 4f; // 1.1 //[1/T

	// Threshold public good concentration before benefit is received
	private static float goodThreshold = 2e-3f;

	// Rate of goods production, to fix by ESS analysis
	private static float goodRateCost1 = 0f; // [1/T]

	private static float goodRateCost2 = 0f; // [1/T]

	private static float costScale = 2f; // 1e-5f; // scales the cost of the

	// public good

	// Yield of biomass on substrate
	private static float Yxs = 0.5f; // [gCOD-PHB/gCOD-S]

	// //////////////////////////////////////////////////////////////////////////////////
	// Solute and Biomass Parameters

	// Public Good
	protected static float goodBulkConcentration = 0f; // [gGood/L]

	// protected static float goodDiffusivity = 4; //4e5f; // [um2/h]
	protected static float goodDiffusivity = 4e5f; // [um2/h]

	// Substrate (Bulk concentration set by input)
	private static float substrateDiffusivity = 4e4f; // [um2/h]

	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	// //////////////////////////////////////////////////////////////////////////////////
	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 250; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = (float) 25 / systemSize;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber1 = 0;

	protected static int initialParticleNumber2 = 0;

	// protected static float maxBiovolume = 50000f; // [um]

	// outpute (write results to file) every:
	// protected static float outputEvery; // .1 [h]

	// /END OF PARAMETERS

	/**
	 * Define bacterial species, solute species, and reaction processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// 2a. Create the solute species (public good only)
		// *******************************************************************************
		// substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));

		// public good
		SoluteSpecies pGood = new SoluteSpecies("pGood", goodDiffusivity);
		pGood.setBulkConcentration(new ConstantBulkConcentration(
				goodBulkConcentration));

		// *******************************************************************************
		// Strain One, quorum sensing
		// Active mass of strain 1, color to be overridden
		ParticulateSpecies activeOne = new ParticulateSpecies("activeOne",
				specificMassBiomass, Color.blue);
		ParticulateSpecies[] spOne = { activeOne };

		float[] fractionalCompositionOne = { 1.0f };

		BiomassSpecies speciesOne = new BiomassSpecies("speciesOne", spOne,
				fractionalCompositionOne);
		speciesOne.setActiveMass(activeOne);

		// ///////////////////////////////////////////////////////////////////////
		// Strain two, constitutive
		// Active mass of strain 2
		ParticulateSpecies activeTwo = new ParticulateSpecies("activeTwo",
				specificMassBiomass, Color.red);

		ParticulateSpecies[] spTwo = { activeTwo };
		float[] fractionalCompositionTwo = { 1.0f };

		BiomassSpecies speciesTwo = new BiomassSpecies("speciesTwo", spTwo,
				fractionalCompositionTwo);
		speciesTwo.setActiveMass(activeTwo);

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// *****************************************************************************

		ProcessFactor substrateUtil = new Saturation(substrate, Ks);
		ProcessFactor goodUtil = new StepReturnSubstrateAndGood(substrate,
				pGood, Ks, goodThreshold);
		ProcessFactor goodUtilGradientsOff = new StepReturnGoodOnly(pGood,
				goodThreshold);

		// 5. Create growth and goods production reactions for each strain
		// *****************************************************************************
		// Strain 1
		// Growth (alternative forms for nutrient gradients ON or OFF
		Reaction growthOne;
		if (gradientsOn) {
			growthOne = new Reaction("growthOne", activeOne, uMax, 1);
			growthOne.addFactor(substrateUtil);
		} else
			growthOne = new Reaction("growthOne", activeOne, uMax, 0);

		// Public Good Production
		Reaction goodProductionOne;
		if (gradientsOn) {
			goodProductionOne = new Reaction("goodProductionOne", activeOne,
					goodRateCost1, 1);
			goodProductionOne.addFactor(substrateUtil);
		} else {
			goodProductionOne = new Reaction("goodProductionOne", activeOne,
					goodRateCost1, 0);
		}

		// Public good use
		Reaction goodUseOne;
		if (gradientsOn) {
			goodUseOne = new Reaction("goodUseOne", activeOne, benefit, 1);
			goodUseOne.addFactor(goodUtil);
		} else {
			goodUseOne = new Reaction("goodUseOne", activeOne, benefit, 1);
			goodUseOne.addFactor(goodUtilGradientsOff);
		}
		// ////////////////////////////////////////////////////////////////////////////

		// Strain 2
		// Growth - alternative forms for nutrient graidents ON or OFF
		Reaction growthTwo;
		if (gradientsOn) {
			growthTwo = new Reaction("growthTwo", activeTwo, uMax, 1);
			growthTwo.addFactor(substrateUtil);
		} else
			growthTwo = new Reaction("growthTwo", activeTwo, uMax, 0);

		// Public good production
		Reaction goodProductionTwo;
		if (gradientsOn) {
			goodProductionTwo = new Reaction("goodProductionTwo", activeTwo,
					goodRateCost2 * uMax, 1);
			goodProductionTwo.addFactor(substrateUtil);
		} else {
			goodProductionTwo = new Reaction("goodProductionTwo", activeTwo,
					goodRateCost2 * uMax, 0);
		}

		// Public Good Use
		Reaction goodUseTwo;
		if (gradientsOn) {
			goodUseTwo = new Reaction("goodUseTwo", activeTwo, benefit, 1);
			goodUseTwo.addFactor(goodUtil);
		} else {
			goodUseTwo = new Reaction("goodUseTwo", activeTwo, benefit, 1);
			goodUseTwo.addFactor(goodUtilGradientsOff);
		}

		// 6. Assign reactions to the species through ReactionStoichiometries
		// ******************************************************************************
		// Strain 1, quorum sensing
		NetReaction rsActiveOne = new NetReaction(3);
		rsActiveOne.addReaction(growthOne, 1);
		rsActiveOne.addReaction(goodUseOne, 1);
		rsActiveOne.addReaction(goodProductionOne, -costScale);

		activeOne.setProcesses(rsActiveOne);

		// Strain 2, constitutive
		NetReaction rsActiveTwo = new NetReaction(3);
		rsActiveTwo.addReaction(growthTwo, 1);
		rsActiveTwo.addReaction(goodUseTwo, 1);
		rsActiveTwo.addReaction(goodProductionTwo, -costScale);

		activeTwo.setProcesses(rsActiveTwo);

		// assign reaction stoichiometry to the solute(s)
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(growthOne, -1 / Yxs);
		rsSubstrate.addReaction(growthTwo, -1 / Yxs);
		substrate.setProcesses(rsSubstrate);

		// public good (G)
		NetReaction rsPgood = new NetReaction(2);
		rsPgood.addReaction(goodProductionOne, 1);
		rsPgood.addReaction(goodProductionTwo, 1);
		pGood.setProcesses(rsPgood);

		// 7. add solute species and biomass species to the system
		addBiomassSpecies(speciesOne);
		addBiomassSpecies(speciesTwo);

		addSoluteSpecies(substrate);
		addSoluteSpecies(pGood);

	}

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumber1, initialParticleNumber2 };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(kdetach);
		setDetachmentHandler(df);
	}

	/**
	 * @param app
	 */
	protected static void setSystemParametersAndInitializeSystemSpace(
			ApplicationComponent app) {
		// create the space
		app.setSystemSpaceParameters(geometry, systemSize,
				relativeMaximumRadius, relativeMinimumRadius,
				relativeBoundaryLayer, gridSide, kShov);
		// initialize
		app.initializeSystemSpace();
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line
		if (args.length < 5) {
			throw new RuntimeException(
					"input arguments missing: \n"
							+ "1: output directory (CAUTION!!! directory will be erased \n"
							+ "2: seed for random number generator \n"
							+ "3: flag for running with graphics (1 on, 0 off) \n"
							+ "4: number of cells at inoculum \n"
							+ "5: flag for nutrient condition (1:lowSubstrate, 2:highSubstrate, 3: noDiff");
		}
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);

		// graphics on/off
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);

		// public goods turned off
		goodRateCost1 = 0;
		goodRateCost2 = 0;

		initialParticleNumber1 = Integer.parseInt(args[3]);
		initialParticleNumber2 = initialParticleNumber1;

		// parse the nutrient concentration
		// 2e-3f produces towers
		// 2e-2f produces sectors?
		// 2e-1f produces nada
		substrateBulkConcentration = Float.parseFloat(args[4]); // [gO2/L] DO 20
		gradientsOn = true;
		// grtadrients off with input substrate 0
		if (substrateBulkConcentration == 0) {
			substrateBulkConcentration = 2e-1f;
			gradientsOn = false;
		}
		// whatever is reached first will stop the simulation
		simulationFinishTime = 500; // [h]
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new ESSGoodsCleanRadial3D();
		// the biovolume series
		VariableSeries biovolume = new BiovolumeSeries();
		VariableSeries[] runLengthSeries = { new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries() };
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// bulk concentrations
			app = new BulkConcentrationVizualizer(app);
			// finally, the controller must be the last decorator to add
			app = new VizualModelControler(app);
		}
		try {
			// create the space
			setSystemParametersAndInitializeSystemSpace(app);
			// initialize
			app.intializeStateWriters(outputDirectory);
			// Pov witer is added twice
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayShovingWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new SoluteConcentrationWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new SolidsConcentrationWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new ParticlePositionWriter()));
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writer
			SimulationResultsWriter writer = new SimulationResultsWriter();
			writer.addSeries(biovolume);
			// spw.addSeries(Model.model().detachedBiomassContainer()
			// .getTotalDetachedBiomassSeries());
			app.addStateWriter(writer);
			// add the time constraints writer
			app.addStateWriter(new TimeConstraintsWriter());
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// add bulk concentrations of all solutes as variable series
			for (Iterator i = Model.model().getSoluteSpecies().iterator(); i
					.hasNext();) {
				SoluteSpecies s = (SoluteSpecies) i.next();
				writer.addSeries(s.getBulkConcentrationSeries());
				writer.addSeries(s.getRateTimeSeries());
			}
			// add particulate global masses to write
			for (Iterator i = Model.model().getParticulateSpecies().iterator(); i
					.hasNext();) {
				ParticulateSpecies s = (ParticulateSpecies) i.next();
				writer.addSeries(s.getTotalMassSeries());
			}
		} catch (ModelException e) {
			System.out.println(e);
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// set writing to infinite to prevent any writing during sim
			//Model.model().setCompulsoryTimeStep(Float.POSITIVE_INFINITY);
			Model.model().setFinishIterationTime(simulationFinishTime);
			Model.model().setMaxBiovolume(2e5f);
			// start the iteration
			app.writeState(); // write iteration 0
			app.startIterating();
		} catch (Exception e1) {
			try {
				app.forceWriteState();
			} catch (IOException e2) {
				System.err.println("Error serializing state:");
				System.err.println("");
				e2.printStackTrace();
			}
			System.err.println("");
			System.err.println("Program failed due to :");
			System.err.println("");
			e1.printStackTrace();
			System.out.println(e1);
		}
		try {
			app.forceWriteState();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}
