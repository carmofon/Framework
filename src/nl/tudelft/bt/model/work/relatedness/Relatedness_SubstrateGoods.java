package nl.tudelft.bt.model.work.relatedness;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.GranuleBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.NoBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Simulates growth of two strains with neutral advantage in a 2D colony
 * spreading on agar plate
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class Relatedness_SubstrateGoods extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	// output directory name
	protected static String outputDirectory = "/Users/careynadell/results/drift2D/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	protected static boolean gradientsOn;
	
	// initial number of particles in the system (inoculum)
	protected static int initialCheater = 25;
	protected static int initialCooperator = 25;
	
	// 1. Reaction parameters
	// *********************************************************************************

	// Effect of public goods upon growth rate
	private static float beta = 10f; // [dim/less]
	
	private static float maintenanceFraction = 0.017f;
	
	private static float basalGrowth = 0.5f;

	// Rate of goods production
	private static float gamma = 1f; // [1/T]

	// cost of producing public goods
	private static float Cg = .1f; // [Mg/Mx]

	// Yield of biomass on substrate
	private static float Yxs = 0.5f; // [gCOD-PHB/gCOD-S]

	// //////////////////////////////////////////////////////////////////////////////////
	// Public Good and Biomass Parameters
	protected static float goodBulkConcentration = 0; // [gGood/L]

	private static float goodDiffusivity = 4e2f; // 4e5f [um2/h]

	protected static float substrateBulkConcentration = 3e-1f; // [gO2/L] DO 20

	private static float substrateDiffusivity = 4e5f; // [um2/h]

	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	// //////////////////////////////////////////////////////////////////////////////////
	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 1500; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 2f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = (float) 50 / systemSize;

	protected static float maximumColonyRadius = 800; // [um]

	protected static float initialColonyRadius = 5; // [um]

	// other model parameters
	protected static int gridSide = 129; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]
	
	// iteration finish time
	protected static float simulationFinishTime = 100f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 0.1f; // [h]

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
		ParticulateSpecies activeCooperator = new ParticulateSpecies("activeCooperator",
				specificMassBiomass, Color.green);
		// Active mass of red strain
		ParticulateSpecies activeCheater = new ParticulateSpecies("activeCheater",
				specificMassBiomass, Color.red);

		// 3. Define particulate compsition of all strains (no polymer)
		// *******************************************************************************
		// green strain
		ParticulateSpecies[] spCooperator = { activeCooperator };
		float[] fractionalVolumeCompositionCooperator = { 1.0f };

		BiomassSpecies speciesCooperator = new BiomassSpecies("speciesCooperator",
				spCooperator, fractionalVolumeCompositionCooperator);
		speciesCooperator.setActiveMass(activeCooperator);

		// red strain
		ParticulateSpecies[] spCheater = { activeCheater };
		float[] fractionalVolumeCompositionCheater = { 1.0f };

		BiomassSpecies speciesCheater = new BiomassSpecies("speciesCheater", spCheater,
				fractionalVolumeCompositionCheater);
		speciesCheater.setActiveMass(activeCheater);


		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// *****************************************************************************

		ProcessFactor gUtil = new Linear2(substrate, pGood);
		ProcessFactor gUtilWithMaintenance = new Linear2MaintBasalGrowth(substrate, pGood,
												maintenanceFraction, basalGrowth);

		// 5. Create growth and goods production reactions for each strain
		// *****************************************************************************
		// Cooperator
		Reaction growthCooperator;
		if (gradientsOn) {
			growthCooperator = new Reaction("growthCooperator", activeCooperator, beta, 1);
			growthCooperator.addFactor(gUtilWithMaintenance);
		}
		else
			growthCooperator = new Reaction("growthCooperator", activeCooperator, beta, 0);
		
		Reaction nutrientUptakeCooperator = new Reaction("nutrientUptakeCooperator", activeCooperator, beta, 1);
		nutrientUptakeCooperator.addFactor(gUtil);

		Reaction pGoodProduction_Cooperator = new Reaction("pGoodProduction_Cooperator",
				activeCooperator, gamma, 0);
		///////////////////////////////////////////////////////////////////////////////
		// Cheater
		Reaction growthCheater;
		if (gradientsOn) {
			growthCheater = new Reaction("growthCheater", activeCheater, beta, 1);
			growthCheater.addFactor(gUtilWithMaintenance);		
		}
		else
			growthCheater = new Reaction("growthCheater", activeCheater, beta, 0);
			
		Reaction nutrientUptakeCheater = new Reaction("nutrientUptakeCheater", activeCheater, beta, 1);
		nutrientUptakeCheater.addFactor(gUtil);
		//////////////////////////////////////////////////////////////////////////////
		
		// Flux of public good
		Reaction fluxPg = new Flux("fluxPg", pGood, goodBulkConcentration);

		// Flux of substrate
		Reaction fluxSubstrate = new Flux("fluxSubstrate", substrate,
				substrateBulkConcentration);

		// 6. Assign reactions to the species through ReactionStoichiometries
		// ******************************************************************************
		// Cooperator
		NetReaction rsActiveGreen = new NetReaction(2);
		rsActiveGreen.addReaction(growthCooperator, 1);
		rsActiveGreen.addReaction(pGoodProduction_Cooperator, -Cg);
		activeCooperator.setProcesses(rsActiveGreen);

		// Cheater
		NetReaction rsActiveRed = new NetReaction(1);
		rsActiveRed.addReaction(growthCheater, 1);
		activeCheater.setProcesses(rsActiveRed);

		// assign reaction stoichiometry to the solute(s)
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(3);
		rsSubstrate.addReaction(nutrientUptakeCooperator, -1 / Yxs);
		rsSubstrate.addReaction(nutrientUptakeCheater, -1 / Yxs);
		
		rsSubstrate.addReaction(fluxSubstrate, substrateDiffusivity
				/ ExtraMath.sq(relativeBoundaryLayer * systemSize));
		substrate.setProcesses(rsSubstrate);

		// public good (G)
		NetReaction rsPgood = new NetReaction(2);
		rsPgood.addReaction(pGoodProduction_Cooperator, 1);
		
		rsPgood.addReaction(fluxPg, goodDiffusivity
				/ ExtraMath.sq(relativeBoundaryLayer * systemSize));
		pGood.setProcesses(rsPgood);

		// 7. add solute species and biomass species to the system
		addBiomassSpecies(speciesCooperator);
		addBiomassSpecies(speciesCheater);
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
		int[] nCells = { initialCooperator, initialCheater };
		inoculateRandomlyInsideRadius(nCells, initialColonyRadius);
	}
	
	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		_boundaryLayer = new NoBoundaryLayer();
		// create the boundary conditions
		MultigridVariable
				.setBoundaryConditions(new GranuleBoundaryConditions());
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(kdetach);
		setDetachmentHandler(df);
		// set the maximum granule radius
		try {
			turnShedingOn();
		} catch (InvalidValueException e) {
			System.out.println(e);
			System.exit(-1);
		}
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
		if (args.length < 4) {
			throw new RuntimeException("input arguments missing:\n"
					+ "1: output directory (CAUTION!!!"
					+ " directory will be erased\n"
					+ "2: seed for random number generator\n"
					+ "3: flag for running with graphics (1 on, 0 off)"
					+ "4: flag for running with nutrient gradients on or off"
					+ "5: number of initial cheater cells [0,50]");
		}
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		gradientsOn = (Integer.parseInt(args[3]) == 1);
		initialCheater = Integer.parseInt(args[4]);
		initialCooperator = 50 - initialCheater;
		System.out.println(initialCheater);
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Relatedness_SubstrateGoods();
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
			spw.addSeries(runLengthSeries[0]);
			spw.addSeries(runLengthSeries[1]);
			spw.addSeries(runLengthSeries[2]);
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