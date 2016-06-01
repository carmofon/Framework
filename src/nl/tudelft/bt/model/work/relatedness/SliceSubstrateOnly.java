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

/**
 * Simulates growth of two strains with neutral advantage in a 2D colony
 * spreading on agar plate
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class SliceSubstrateOnly extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program


	// output directory name
	protected static String outputDirectory = "/Users/careynadell/results/drift2D/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	protected static boolean gradientsOn;
	
	// 1. Reaction parameters
	// *********************************************************************************

	//Maintenance Fraction
	private static float maintenanceFraction = 0.015f;
	
	//max growth rate on substrate alone
	private static float uMax = 1.5f;
	
	// Effect of public goods upon growth rate
	private static float benefit = 1.1f; // [1/T

	// Rate of goods production, to fix by ESS analysis
	private static float goodRate = 1f; // [1/T]

	// cost of producing public goods
	private static float cost = 0.41f; // [Mg/Mx]

	// Yield of biomass on substrate
	private static float Yxs = 0.5f; // [gCOD-PHB/gCOD-S]

	// //////////////////////////////////////////////////////////////////////////////////
	// Public Good and Biomass Parameters
	protected static float goodBulkConcentration = 0; // [gGood/L]

	private static float goodDiffusivity = 4e3f; // 4e5f [um2/h]

	protected static float substrateBulkConcentration = 8e-2f; //2e-1f; // [gO2/L] DO 20

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
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = (float) 25 / systemSize;

	//protected static float maximumColonyRadius = 800; // [um]

	//protected static float initialColonyRadius = 5; // [um]

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber1 = 100;
	protected static int initialParticleNumber2 = 100;
	
	// iteration finish time
	protected static float simulationFinishTime = 400f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 1f; // [h]

	// /END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// 2a. Create the solute species (public good only)
		// *******************************************************************************
		// Public Good (G)
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));

		SoluteSpecies pGood = new SoluteSpecies("pGood", goodDiffusivity);
		pGood.setBulkConcentration(new ConstantBulkConcentration(
				goodBulkConcentration));
		//

		// 2b. Create the particulate species (strains)
		// ******************************************************************************
		// Active mass of green strain
		ParticulateSpecies activeGreen = new ParticulateSpecies("activeGreen",
				specificMassBiomass, Color.green);
		// Active mass of red strain
		ParticulateSpecies activeRed = new ParticulateSpecies("activeRed",
				specificMassBiomass, Color.red);

		// 3. Define particulate compsition of all strains (active mass only, no
		// polymer)
		// *******************************************************************************
		// green strain
		ParticulateSpecies[] spGreen = { activeGreen };
		float[] fractionalVolumeCompositionGreen = { 1.0f };

		BiomassSpecies speciesGreen = new BiomassSpecies("speciesGreen",
				spGreen, fractionalVolumeCompositionGreen);
		speciesGreen.setActiveMass(activeGreen);

		// red strain
		ParticulateSpecies[] spRed = { activeRed };
		float[] fractionalVolumeCompositionRed = { 1.0f };

		BiomassSpecies speciesRed = new BiomassSpecies("speciesRed", spRed,
				fractionalVolumeCompositionRed);
		speciesRed.setActiveMass(activeRed);


		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// *****************************************************************************

		ProcessFactor gUtil = new Linear1(substrate);
		ProcessFactor gUtilMaintenance = new Linear1WithMaintenance(substrate, maintenanceFraction);

		// 5. Create growth and goods production reactions for each strain
		// *****************************************************************************
		// Green
		Reaction growthGreen;
		if (gradientsOn) {
			growthGreen = new Reaction("growthGreen", activeGreen, uMax, 1);
			growthGreen.addFactor(gUtilMaintenance);
		}
		else
			growthGreen = new Reaction("growthGreen", activeGreen, uMax, 0);
		
		Reaction nutrientUptakeGreen = new Reaction("nutrientUptakeGreen", activeGreen, uMax, 1);
		nutrientUptakeGreen.addFactor(gUtil);

		//Reaction pGoodProduction_Green = new Reaction("pGoodProduction_Green", activeGreen, 0, 0);

		// Red
		Reaction growthRed;
		if (gradientsOn) {
			growthRed = new Reaction("growthRed", activeRed, uMax, 1);
			growthRed.addFactor(gUtilMaintenance);		
		}
		else
			growthRed = new Reaction("growthRed", activeRed, uMax, 0);
			
		Reaction nutrientUptakeRed = new Reaction("nutrientUptakeRed", activeRed, uMax, 1);
		nutrientUptakeRed.addFactor(gUtil);
		
		//Reaction pGoodProduction_Red = new Reaction("pGoodProduction_Red", activeRed, 0, 0);

		// 6. Assign reactions to the species through ReactionStoichiometries
		// ******************************************************************************
		// active Green
		NetReaction rsActiveGreen = new NetReaction(1);
		rsActiveGreen.addReaction(growthGreen, 1);
		//rsActiveGreen.addReaction(pGoodProduction_Green, -cost);
		activeGreen.setProcesses(rsActiveGreen);

		// active Red
		NetReaction rsActiveRed = new NetReaction(1);
		rsActiveRed.addReaction(growthRed, 1);
		//rsActiveRed.addReaction(pGoodProduction_Red, -cost);
		activeRed.setProcesses(rsActiveRed);

		// assign reaction stoichiometry to the solute(s)
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(nutrientUptakeGreen, -1 / Yxs);
		rsSubstrate.addReaction(nutrientUptakeRed, -1 / Yxs);
	
		substrate.setProcesses(rsSubstrate);

		// public good (G)
		//NetReaction rsPgood = new NetReaction(6);
		//rsPgood.addReaction(pGoodProduction_Green, 1);
		//rsPgood.addReaction(pGoodProduction_Red, 1);
		
		//pGood.setProcesses(rsPgood);

		// 7. add solute species and biomass species to the system
		addBiomassSpecies(speciesGreen);
		addBiomassSpecies(speciesRed);
		
		addSoluteSpecies(substrate);
		//addSoluteSpecies(pGood);

	}

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumber1, initialParticleNumber2};
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
		if (args.length < 3) {
			throw new RuntimeException("input arguments missing:\n"
					+ "1: output directory (CAUTION!!!"
					+ " directory will be erased\n"
					+ "2: seed for random number generator\n"
					+ "3: flag for running with graphics (1 on, 0 off)"
					+ "4: flag for running with nutrient gradients on or off");
		}
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		gradientsOn = (Integer.parseInt(args[3]) == 1);
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new SliceSubstrateOnly();
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
			Model.model().setCompulsoryTimeStep(outputEvery);
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