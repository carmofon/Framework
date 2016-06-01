package nl.tudelft.bt.model.work.relatedness;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.util.ExtraMath;
import nl.tudelft.bt.model.work.quorumsensing.QSBiomassSpecies;
import nl.tudelft.bt.model.work.quorumsensing.UpRegulationStep;

/**
 * Simulates growth of two strains with neutral advantage in a 2D colony
 * spreading on agar plate
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class QS_ESSgoods extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/careynadell/results/QSpractice/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	// These are all now defined in the main function
	protected static boolean gradientsOn;

	protected static float simulationFinishTime; // [h]

	protected static float substrateBulkConcentration;

	protected static int nutrientCondition;

	// 1. Reaction parameters
	// *********************************************************************************

	// max growth rate on substrate alone
	// private static float uMax = 1.5f;
	private static float uMax = .25f;  //1.25f for samples

	private static float Ks = 3.5e-5f;

	private static float Kg = .1f; // ??

	// Effect of public goods upon growth rate
	private static float benefit; // 1.1 //[1/T
	
	private static float costScale = 0.25f; // scales the cost of the public good

	// Rate of goods production, to fix by ESS analysis
	private static float goodRateCost1 = 0f; // [1/T]

	private static float goodRateCost2 = 0f; // [1/T]	
	
	// Yield of biomass on substrate
	private static float Yxs = 0.5f; // [gCOD-PHB/gCOD-S]
	
	// Yield of Autoinducer on biomass
	private static float Yax = 200f / 2.65e-3f; // [A/gX]

	// //////////////////////////////////////////////////////////////////////////////////
	// Solute and Biomass Parameters
	
	// Public Good
	protected static float goodBulkConcentration = 0f; // [gGood/L]

	protected static float goodDiffusivity = 4e3f; // 4e5f [um2/h]
	
	// Autoinducer (A)
	protected static float autoinducerBulkConcentration = 0.0f;

	protected static float autoinducerDiffusivity = 1.0e4f;
	
	// Autoinducer production constant
	protected static float aiProductionRate = 0.2f; // [1/h]

	// Autoinducer quorum concentration threshold
	protected static float QSthreshold = 1f; // [1/L^2]


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
	protected static float relativeMaximumRadius = 1f / systemSize;

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
	protected static int initialParticleNumber1;

	protected static int initialParticleNumber2;

	//protected static float maxBiovolume = 40000f; // [um]

	// outpute (write results to file) every:
	protected static float outputEvery; // .1 [h]

	// /END OF PARAMETERS

	/**
	 * Define bacterial species, solute species, and reaction processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// 2a. Create the solute species (public good only)
		// *******************************************************************************
		// substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate", substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(substrateBulkConcentration));
		
		// public good
		SoluteSpecies pGood = new SoluteSpecies("pGood", goodDiffusivity);
		pGood.setBulkConcentration(new ConstantBulkConcentration(goodBulkConcentration));
		
		// autoinducer
		SoluteSpecies autoinducer = new SoluteSpecies("auotinducer", autoinducerDiffusivity);
		autoinducer.setBulkConcentration(new ConstantBulkConcentration(autoinducerBulkConcentration));

		// *******************************************************************************
		// Strain One, quorum sensing
		// Active mass of strain 1, color to be overridden
		ParticulateSpecies activeOne = new ParticulateSpecies(
				"activeOne", specificMassBiomass, Color.green);
		ParticulateSpecies[] spOne = { activeOne };
		
		float[] fractionalCompositionOne = { 1.0f };

		QSBiomassSpecies speciesOne = new QSBiomassSpecies("speciesOne",
				spOne, fractionalCompositionOne, autoinducer,
				QSthreshold, Color.green, Color.red);
		speciesOne.setActiveMass(activeOne);
		
		/////////////////////////////////////////////////////////////////////////
		// Strain two, constitutive
		// Active mass of strain 2
		ParticulateSpecies activeTwo = new ParticulateSpecies(
				"activeTwo", specificMassBiomass, Color.blue);
		
		ParticulateSpecies[] spTwo = { activeTwo };
		float[] fractionalCompositionTwo = { 1.0f };

		BiomassSpecies speciesTwo = new BiomassSpecies("speciesTwo",
				spTwo, fractionalCompositionTwo);
		speciesTwo.setActiveMass(activeTwo);


		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// *****************************************************************************

		ProcessFactor substrateUtil = new Saturation(substrate, Ks);
		ProcessFactor goodUtil = new MonodTwoSubstrates(substrate, Ks, pGood,
				Kg);
		ProcessFactor goodUtilGradientsOff = new Saturation(pGood, Kg);
		
		ProcessFactor quorumSensing = new UpRegulationStep(autoinducer, QSthreshold);

		// 5. Create growth and goods production reactions for each strain
		// *****************************************************************************
		// Strain 1
		// Growth (alternative forms for nutrient gradients ON or OFF
		Reaction growthOne;
		if (gradientsOn) {
			growthOne = new Reaction("growthOne",
					activeOne, uMax, 1);
			growthOne.addFactor(substrateUtil);
		} else
			growthOne = new Reaction("growthOne",
					activeOne, uMax, 0);

		// Autoinducer secretion
		Reaction aiProductionOne = new Reaction("aiProductionOne",
				activeOne, aiProductionRate, 0);
		
		// Public Good Production
		Reaction goodProductionOne;
		if (gradientsOn) {
			goodProductionOne = new Reaction("goodProductionOne", activeOne,
					goodRateCost1, 2);
			goodProductionOne.addFactor(substrateUtil);
			goodProductionOne.addFactor(quorumSensing);   // good production under positive QS regulation
		} else {
			goodProductionOne = new Reaction("goodProductionOne", activeOne,
					goodRateCost1, 1);
			goodProductionOne.addFactor(quorumSensing);   // good production under positive QS regulation
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

		// Autoinducer secretion
		Reaction aiProductionTwo = new Reaction("aiProductionTwo",
				activeTwo, aiProductionRate, 0);
		
		// Public good production
		Reaction goodProductionTwo;
		if (gradientsOn) {
			goodProductionTwo = new Reaction("goodProductionTwo", activeTwo,
					goodRateCost2, 1);
			goodProductionTwo.addFactor(substrateUtil);
		} else {
			goodProductionTwo = new Reaction("goodProductionTwo", activeTwo,
					goodRateCost2, 0);
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
		NetReaction rsActiveOne = new NetReaction(4);
		rsActiveOne.addReaction(growthOne, 1);
		rsActiveOne.addReaction(goodUseOne, 1);
		rsActiveOne.addReaction(aiProductionOne, -(1/Yax));
		rsActiveOne.addReaction(goodProductionOne, -costScale);

		activeOne.setProcesses(rsActiveOne);

		// Strain 2, constitutive
		NetReaction rsActiveTwo = new NetReaction(4);
		rsActiveTwo.addReaction(growthTwo, 1);
		rsActiveTwo.addReaction(goodUseTwo, 1);
		rsActiveTwo.addReaction(aiProductionTwo, -(1/Yax));
		rsActiveTwo.addReaction(goodProductionTwo, -costScale);

		activeTwo.setProcesses(rsActiveTwo);

		// assign reaction stoichiometry to the solute(s)
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(growthOne, -1 / Yxs);
		rsSubstrate.addReaction(growthTwo, -1 / Yxs);
		substrate.setProcesses(rsSubstrate);

		// Autoinducer - produced by both strains
		NetReaction rsAutoinducer = new NetReaction(2);
		rsAutoinducer.addReaction(aiProductionOne, 1);
		rsAutoinducer.addReaction(aiProductionTwo, 1);
		autoinducer.setProcesses(rsAutoinducer);
		
		// public good (G)
		NetReaction rsPgood = new NetReaction(4);
		rsPgood.addReaction(goodProductionOne, 1);
		rsPgood.addReaction(goodProductionTwo, 1);
		rsPgood.addReaction(goodUseOne, -1);
		rsPgood.addReaction(goodUseTwo, -1);
		pGood.setProcesses(rsPgood);

		// 7. add solute species and biomass species to the system
		addBiomassSpecies(speciesOne);
		addBiomassSpecies(speciesTwo);

		addSoluteSpecies(substrate);
		addSoluteSpecies(autoinducer);
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
		if (args.length < 6) {
			throw new RuntimeException(
					"input arguments missing: \n"
							+ "1: output directory (CAUTION!!! directory will be erased \n"
							+ "2: seed for random number generator \n"
							+ "3: flag for running with graphics (1 on, 0 off) \n"
							+ "4: good production rate 2 \n"
							+ "5: flag for focal biofilm (=1) or resident biofilm (=0) \n"
							+ "6: flag for nutrient condition (1:lowSubstrate, 2:highSubstrate, 3: noDiff"
							+ "7: benefit of receiving public goods");
		}
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);

		goodRateCost2 = Float.parseFloat(args[3]);
		goodRateCost1 = goodRateCost2 + 0.5f;

		if (Integer.parseInt(args[4]) == 1) {
			initialParticleNumber1 = 12;
			initialParticleNumber2 = 120 - initialParticleNumber1;
		} else {
			initialParticleNumber1 = 120;
			initialParticleNumber2 = 0;
		}
		nutrientCondition = Integer.parseInt(args[5]);

		if (nutrientCondition == 1) {
			// lowSubstrate: Towers
			substrateBulkConcentration = 2.5e-1f; //.5e-1f; // [gO2/L] DO 20
			simulationFinishTime = 200f; // 40 (sample) // [h]
			outputEvery = simulationFinishTime - .01f;
			gradientsOn = true;
		} else {
			// highSubstrate: Sectors
			if (nutrientCondition == 2) {
				substrateBulkConcentration = 1e0f;//0.22e0f; // [gO2/L] DO 20
				simulationFinishTime = 100f; // 15 (Sample) // [h]
				outputEvery = simulationFinishTime - .01f;
				gradientsOn = true;
			} else {
				// No nutrient gradients
				substrateBulkConcentration = 0.5e-1f; // [gO2/L] DO 20
				simulationFinishTime = 10f; // 3.5 (sample) // [h]
				outputEvery = simulationFinishTime - .01f;
				gradientsOn = false;
			}
		}
		
		benefit = Float.parseFloat(args[6]);
		
		float alpha_tilde = 0.0000000000001f;
		float L2 = 100 * 100;
		QSthreshold = alpha_tilde * L2 * aiProductionRate * specificMassBiomass
				/ autoinducerDiffusivity;

		// set numerics for multigrid
		MultigridVariable.setSteps(10, 100);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new QS_ESSgoods();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
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
							new PovRayWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
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
				spw.addSeries(s.getBulkConcentrationSeries());
				spw.addSeries(s.getRateTimeSeries());
			}
			// add particulate global masses to write
			for (Iterator i = Model.model().getParticulateSpecies().iterator(); i
					.hasNext();) {
				ParticulateSpecies s = (ParticulateSpecies) i.next();
				spw.addSeries(s.getTotalMassSeries());
			}
		} catch (ModelException e) {
			System.out.println(e);
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			// Model.model().setCompulsoryTimeStep(outputEvery);
			//Model.model().setMaxBiovolume(maxBiovolume);
			 Model.model().setFinishIterationTime(simulationFinishTime);
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
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}