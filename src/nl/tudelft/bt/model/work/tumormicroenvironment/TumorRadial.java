package nl.tudelft.bt.model.work.tumormicroenvironment;

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
 * Sumlated the gorwth of a tumor with cancer cells and macrophages
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class TumorRadial extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/careynadell/practice";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	// protected static float simulationFinishTime; // [h]

	protected static float finishRadiusRatio = 4.0f; // [um]

	protected static float substrateBulkConcentration;

	protected static float nutrientCondition;

	// 1. Reaction parameters
	// *********************************************************************************

	// max growth rate on substrate alone
	// private static float uMax = 1.5f;
	private static float uMax = 1f;

	private static float Ks = 3.5e-5f;

	private static float Kp = 5e-3f;

	private static float Ki = 10f;

	// Yield of biomass on substrate
	private static float Yxs = 0.5f; // [gCOD-PHB/gCOD-S]

	// //////////////////////////////////////////////////////////////////////////////////
	// Solute and Biomass Parameters

	// Substrate (Bulk concentration set by input)
	private static float substrateDiffusivity = 4e4f; // [um2/h]

	// Particulate species (biomass H)
	protected static float specificMassBiomass = 428f; // [gCOD-H/L]

	// //////////////////////////////////////////////////////////////////////////////////
	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 2000; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 7.5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = (float) 200f / systemSize;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumberCancer = 0;
	protected static int initialParticleNumberMacrophage = 0;
	protected static int initialParticleNumberTCell = 0;

	// protected static float maxBiovolume = 50000f; // [um]

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
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		// IL-4
		SoluteSpecies il4 = new SoluteSpecies("IL-4", substrateDiffusivity);
		il4.setBulkConcentration(new ConstantBulkConcentration(0));
		// Protease
		SoluteSpecies protease = new SoluteSpecies("protease",
				substrateDiffusivity);
		protease.setBulkConcentration(new ConstantBulkConcentration(0));

		// *******************************************************************************
		// Cell type one, cancer cell
		ParticulateSpecies activeCancer = new ParticulateSpecies(
				"activeCancer", specificMassBiomass, Color.red);
		ParticulateSpecies[] spOne = { activeCancer };

		float[] fractionalCompositionOne = { 1.0f };

		BiomassSpecies speciesCancer = new BiomassSpecies("cancer", spOne,
				fractionalCompositionOne);
		speciesCancer.setInducibleColor(protease, Kp, Color.red, new Color(
				0.1f, 0f, 0f));

		// ///////////////////////////////////////////////////////////////////////
		// Cell type two, macrophage
		ParticulateSpecies uninducedMacrophage = new ParticulateSpecies(
				"uninducedMacrophage", specificMassBiomass, Color.gray);

		ParticulateSpecies[] spTwo = { uninducedMacrophage };
		float[] fractionalCompositionTwo = { 1.0f };

		BiomassSpecies speciesMacrophage = new BiomassSpecies("macrophage",
				spTwo, fractionalCompositionTwo);
		// TODO delete
		// speciesMacrophage.setActiveMass(uninducedMacrophage);
		speciesMacrophage.setInducibleColor(il4, Ki, Color.green, new Color(0f,
				0.1f, 0f));

		// ///////////////////////////////////////////////////////////////////////
		// Cell type three, T-cell
		ParticulateSpecies activeTCell = new ParticulateSpecies("activeTCell",
				specificMassBiomass, Color.blue);

		ParticulateSpecies[] spThree = { activeTCell };
		float[] fractionalCompositionThree = { 1.0f };

		BiomassSpecies speciesTCell = new BiomassSpecies("tCell", spThree,
				fractionalCompositionThree);
		speciesTCell.setActiveMass(activeTCell);

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// *****************************************************************************

		ProcessFactor substrateMonod = new Saturation(substrate, Ks);
		ProcessFactor il4Switch = new Saturation(il4, Ki);
		ProcessFactor proteaseSaturation = new Saturation(protease, Kp);

		// 5. Create growth and goods production reactions for each strain
		// *****************************************************************************
		// Cancer growth
		Reaction growthCancerBasal = new Reaction("growthCancerBasal",
				activeCancer, uMax * 0.001f, 0);
		// Cancer IL-4 secretion
		Reaction il4SecretionCancer = new Reaction("il4SecretionCancer",
				activeCancer, uMax * 0.02f, 0);
		//il4SecretionCancer.addFactor(substrateMonod);
		// ////////////////////////////////////////////////////////////////////////////

		// Macrophage growth
		Reaction macrophageProliferation = new Reaction(
				"macrophageProliferation", uninducedMacrophage, uMax*0.5f, 1);
		macrophageProliferation.addFactor(substrateMonod);
		// Protease secretion
		Reaction proteaseSecretion = new Reaction("proteaseSecretion",
				uninducedMacrophage, uMax, 2);
		proteaseSecretion.addFactor(il4Switch);
		proteaseSecretion.addFactor(substrateMonod);

		// T-Cell growth
		Reaction iL4SecretionTCell = new Reaction("tCellIL4Secretion", activeTCell,
				uMax * 0.02f * 4, 0);

		// 6. Assign reactions to the species through ReactionStoichiometries
		// ******************************************************************************
		// cancer
		NetReaction rsActiveCancer = new NetReaction(1);
		rsActiveCancer.addReaction(growthCancerBasal, 1);
		//rsActiveCancer.addReaction(invasionCancer, 1);
		activeCancer.setProcesses(rsActiveCancer);

		// macrophage
		NetReaction rsActiveMacrophage = new NetReaction(1);
		// rsActiveMacrophage.addReaction(growthMacrophageBasal, 1);
		rsActiveMacrophage.addReaction(macrophageProliferation, 1f);
		uninducedMacrophage.setProcesses(rsActiveMacrophage);

		// T-cell
		NetReaction rsActiveTCell = new NetReaction(1);
		rsActiveTCell.addReaction(iL4SecretionTCell, 0);

		activeTCell.setProcesses(rsActiveTCell);

		// assign reaction stoichiometry to the solute(s)
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(1);
		//rsSubstrate.addReaction(invasionCancer, -1 / Yxs);
		rsSubstrate.addReaction(macrophageProliferation, -1 / Yxs);
		// rsSubstrate.addReaction(growthTCell, -1 / Yxs);
		substrate.setProcesses(rsSubstrate);

		// IL-4
		NetReaction rsIl4 = new NetReaction(2);
		rsIl4.addReaction(iL4SecretionTCell, 1f);
		rsIl4.addReaction(il4SecretionCancer, 1f);
		il4.setProcesses(rsIl4);

		// protease
		NetReaction rsProtease = new NetReaction(1);
		rsProtease.addReaction(proteaseSecretion, 1);
		protease.setProcesses(rsProtease);

		// 7. add solute species and biomass species to the system
		addBiomassSpecies(speciesCancer);
		addBiomassSpecies(speciesMacrophage);
		addBiomassSpecies(speciesTCell);

		addSoluteSpecies(substrate);
		addSoluteSpecies(il4);
		addSoluteSpecies(protease);

	}

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumberCancer,
				initialParticleNumberMacrophage, initialParticleNumberTCell };
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
							+ "4: number of cells at inoculum \n"
							+ "5: flag for nutrient condition (1:lowSubstrate, 2:highSubstrate, 3: noDiff"
							+ "6: 2D/3D (2:2D, 3:3D");
		}
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);

		// initial number of cells
		initialParticleNumberCancer = Math
				.round(Float.parseFloat(args[3]) * 0.9f);
		initialParticleNumberMacrophage = Math
				.round(Float.parseFloat(args[3]) * 0.05f);
		initialParticleNumberTCell = Math
				.round(Float.parseFloat(args[3]) * 0.1f);
		//initialParticleNumberTCell = 0;

		int nutrientCondition = Integer.parseInt(args[4]);
		int dimensionallity = Integer.parseInt(args[5]);
		switch (dimensionallity) {
		case 2:
			geometry = 2;
			break;
		case 3:
			geometry = 3;
			break;
		default:
			System.out.println("input geometery (2D/3D");
			System.exit(-1);
			break;
		}

		if (nutrientCondition == 1) {
			// lowSubstrate: Towers
			substrateBulkConcentration = 2e-1f; // [gO2/L] DO 20
		} else {
			substrateBulkConcentration = .5e0f; // [gO2/L] DO 20
		}

		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TumorRadial();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorating the application
			app = new BiomassVizualizer(app);
		}
		RunLengthXSeries x = new RunLengthXSeries();
		RunLengthXSeries y = new RunLengthXSeries();
		try {
			// create the space
			setSystemParametersAndInitializeSystemSpace(app);
			// initialize
			app.intializeStateWriters(outputDirectory);
			// Pov writer is added twice
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
			// the simulation parameters writer
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(new BiovolumeSeries());
			// set up the run length
			spw.addSeries(x);
			spw.addSeries(y);
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
			// makes sure it never writes outputs during sim
			// Model.model().setCompulsoryTimeStep(1.0f);
			// Model.model().setMaxBiovolume(maxBiovolume);
			// Model.model().setFinishIterationTime(simulationFinishTime);
			// calulate the final radius as 4 x the initial radiius
			float xValue = x.getLastY() * 0.5f;
			float yValue = y.getLastY() * 0.5f;
			float finishSimulationRadius = (xValue * 0.5f + yValue * 0.5f)
					* finishRadiusRatio;
			// Model.model().setMaxRunLength(finishSimulationRadius*2.0f);
			Model.model().setMaxRunLength(1000.0f);
			// start the iteration
			// app.writeState(); // do not write iteration 0
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
			System.out.println("Writing the output at end of simulation");
			app.forceWriteStateWithoutSerializing();
		} catch (ModelException e2) {
			System.err.println("Error trying to write output");
		}
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}
