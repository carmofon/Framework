package nl.tudelft.bt.model.work.carlos;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
//import nl.tudelft.bt.model.apps.ApplicationComponentExtension1;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_conditions.GranuleBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.FixedOnTopBoundaryLayer;
import nl.tudelft.bt.model.multigrid.boundary_layers.NoBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassSpecies;
import nl.tudelft.bt.model.particlebased.tumor.TumorModelHandler;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.work.drift.BiomassWithNeutralColors;

/**
 * We use two populations P1 and P2 which cooperates at different length scales.
 * Both populations grow on a fraction f of the external nutrient S while
 * internalize the rest (1-f)S proportionally to the local concentration of
 * growh factor G. The growth factor G diffuses in the media.
 * 
 * @author Vanni Bucci (bucciv@mskcc.org)
 */
public class ShortLongCooperation extends TumorModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/CarmonaFontaine/Documents/Xlab/Models/Cooperation";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// Nutrient flux parameters
	protected static float externalMassTransferCoefficient = 1e4f; // [h-1]
	// protected static float externalMassTransferCoefficient = 1e-4f; // [h-1]

	// Reaction parameters
	// GROWTH
	protected static float uMax = 1f; // [h-1]
	private static float KS = 10f; // [gS/L] Merle
	private static float YS = 1f; // [gX/gS]
	// fraction of nutrient readily internalizable
	private static float f = 0.5f;

	// GROWTH FACTOR
	// half saturation constant for gfs
	private static float KG = 1f; // [gG/L]
	// gfs yield
	private static float Yg = 1f; // [gG/gX];
	// gfs prod rate
	protected static float qg1 = 10f; // [gX/m3 s]
	protected static float qg2 = 10f; // [gX/m3 s]
	// gf1 diffusivity
	private static float Dg1 = 1f; // [m2/s]
	// gf2 diffusivity
	private static float Dg2 = 1000f; // [m2/s]
	// G1 & G2
	protected static float g1BulkConcentration = 0f; // [gS/L] DO 20
	protected static float g2BulkConcentration = 0f; // [gS/L] DO 20

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	// SUBSTRATE
	protected static float sBulkConcentration = 20f; // [gS/L] DO 20
	private static float substrateDiffusivity = 1f; // [m2/s]

	// SIZE PARAMETERS
	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 25	; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 0 / systemSize;

	// This turns-off cell shedding from colony
	protected static float maximumColonyRadius = Float.POSITIVE_INFINITY; // [um]

	protected static float initialColonyRadius = 1; // [um]

	// other model parameters
	protected static int gridSide = 17; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 200;

	// iteration finish time
	protected static float simulationFinishTime = 100f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 0.15f; // [h]

	private static SoluteSpecies substrate;
	private static SoluteSpecies G1;
	private static SoluteSpecies G2;

	// /END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {

		// ---------------------------------------------------
		// ---------------------------------------------------
		// 1. Create the solute species
		// substrate (S)
		substrate = new SoluteSpecies("substrate", substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				sBulkConcentration));
		// Toxin species
		G1 = new SoluteSpecies("G1", Dg1);
		// set up the simplest type of bulk concentration: constant
		G1.setBulkConcentration(new ConstantBulkConcentration(
				g1BulkConcentration));
		G2 = new SoluteSpecies("G2", Dg2);
		// set up the simplest type of bulk concentration: constant
		G2.setBulkConcentration(new ConstantBulkConcentration(
				g2BulkConcentration));

		// ---------------------------------------------------
		// ---------------------------------------------------
		// 2. Create the particulate species (solids)
		// The active mass of P1 strain
		ParticulateSpecies activeP1 = new ParticulateSpecies("activeP1",
				specificMassBiomass, Color.green);
		// P1 strain
		ParticulateSpecies[] spP1 = { activeP1 };
		float[] fractionalVolumeCompositionP1 = { 1.0f };

		// The active mass of P2 strain
		ParticulateSpecies activeP2 = new ParticulateSpecies("activeP2",
				specificMassBiomass, Color.red);
		// P2 strain
		ParticulateSpecies[] spP2 = { activeP2 };
		float[] fractionalVolumeCompositionP2 = { 1.0f };
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)

		// ---------------------------------------------------
		// ---------------------------------------------------
		// 3. Create the biomass species
		// P1
		BiomassSpecies P1 = new GranuleBiomassSpecies("speciesP1", spP1,
				fractionalVolumeCompositionP1);
		P1.setActiveMass(activeP1);
		//P1.setInducibleColor(G1, KG, Color.green,
		//		new Color(0, 0.5f, 0));
		
		// P2
		BiomassSpecies P2 = new GranuleBiomassSpecies("speciesP2", spP2,
				fractionalVolumeCompositionP2);
		P2.setActiveMass(activeP2);
		//P2.setInducibleColor(G2, KG, Color.red,
		//		new Color(0.5f, 0, 0));

		// ---------------------------------------------------
		// ---------------------------------------------------
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// ProcessFactor mS = new Saturation(oxygen, KO);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		// for XH
		// ProcessFactor mS = new Saturation(substrate, KS);
		// ProcessFactor mg = new Monod_Sum_of_TwoSubstrates(G1,G2, KG);
		// ProcessFactor mg2 = new Saturation(G2, KG);
		ProcessFactor mS = new S_effective_TG(substrate, G1, G2, KS, KG, f);

		// ProcessFactor mSWithMaintenace = new SaturationWithMaintenance(
		// substrate, KS, fMaintenance);

		// ---------------------------------------------------
		// ---------------------------------------------------
		// 5. Create the reactions
		// growth P1
		Reaction GrowthP1 = new Reaction("growthP1", activeP1, uMax, 1);
		GrowthP1.addFactor(mS);
		// aerobic growth P2
		Reaction GrowthP2 = new Reaction("growthP2", activeP2, uMax, 1);
		GrowthP2.addFactor(mS);

		// g1 production
		Reaction G1Production = new Reaction("g1Production", activeP1, qg1, 1);
		// g2 production
		Reaction G2Production = new Reaction("g2Production", activeP2, qg2, 1);
		//
		Reaction flux = new Flux("flux", substrate, sBulkConcentration);
		Reaction fluxg1 = new Flux("flux", G1, g1BulkConcentration);
		Reaction fluxg2 = new Flux("flux", G2, g2BulkConcentration);
		// ---------------------------------------------------
		// ---------------------------------------------------
		// 6. Assign reaction to the species through Reaction Stoichiometries
		// activeP1
		NetReaction rsXP1 = new NetReaction(1);
		rsXP1.addReaction(GrowthP1, 1);
		activeP1.setProcesses(rsXP1);
		// activeP2
		NetReaction rsXP2 = new NetReaction(1);
		rsXP2.addReaction(GrowthP2, 1);
		activeP2.setProcesses(rsXP2);
		// assign reaction stoichiometry to the solutes
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(3);
		rsSubstrate.addReaction(GrowthP1, -1 / YS);
		rsSubstrate.addReaction(GrowthP1, -1 / YS);
		rsSubstrate.addReaction(flux, externalMassTransferCoefficient);
		substrate.setProcesses(rsSubstrate);
		// g1
		NetReaction rsG1 = new NetReaction(2);
		rsG1.addReaction(G1Production, Yg);
		rsG1.addReaction(fluxg1, externalMassTransferCoefficient);
		G1.setProcesses(rsG1);
		// g1
		NetReaction rsG2 = new NetReaction(2);
		rsG2.addReaction(G2Production, Yg);
		rsG2.addReaction(fluxg2, externalMassTransferCoefficient);
		G2.setProcesses(rsG2);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(P1);
		addBiomassSpecies(P2);
		addSoluteSpecies(substrate);
		addSoluteSpecies(G1);
		addSoluteSpecies(G2);
	}

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
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
	protected void inoculate() {
		int[] nCells = { initialParticleNumber, initialParticleNumber };
		// inoculateRandomlyInsideRadius(nCells, initialColonyRadius);
		inoculateRandomlyEverywhere(nCells, initialColonyRadius);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(kdetach);
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
					+ "3: flag for running with graphics (1 on, 0 off)");
		}
		// set external mass transfer, the 4th, optional, input
		if (args.length > 3) {
			externalMassTransferCoefficient = Float.parseFloat(args[3]);
		}

		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new ShortLongCooperation();
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
			// initialize
			app.initializeDiffusionReactionSystem(); // also inoculates
			// Pov writer is added twice
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new PovRayWriter()));
			app.addTimedStateWriter(new ImageWriter(substrate));
			app.addTimedStateWriter(new ImageWriter(G1));
			app.addTimedStateWriter(new ImageWriter(G2));
			app.addTimedStateWriter(new ImageWriter());
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
			//app.startIterating2(true, 0.1f);
			//app.KillByStochasticEffect(0.1f);
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