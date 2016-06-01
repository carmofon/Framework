package nl.tudelft.bt.model.work.epsproducers.mutation;

import java.awt.Color;
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

/**
 * Determine the impact of a single change from EPS+ to EPS- on the neighborhood
 * Keeping distances to the mutant
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class EpsProducersMutationKeepingDistances extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	static String outputDirectory = "c:/results/EpsProducersMutationImpactDistancesNoMut/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "e:/results/example1/" directory "e:\results" must
	// exist in your computer. The subdirectory "example1" will be created
	// if it is non-existant. If it exists, its contents will be deleted
	// during the program initialization

	// geometry (default is 2D) - change value to 3 for 3D
	protected static int geometry = 2;

	// Solute species
	// Substrate (S) - the only solute species used here
	// protected static float oxygenBulkConcentration = 8e-3f * 0.5f; // [g/L]
	protected static float oxygenBulkConcentration = 8e-3f; // [g/L]

	// protected static float oxygenBulkConcentration = 8f; // [g/L]

	// 
	static float oxygenDiffusivity = 8.33e6f; // [um^2/h]

	// 

	//
	// Particulate species
	// biomass X
	// protected static float specificMassX = 93.3f; // [g-C/L] (Bratbak 1985)
	protected static float specificMassX = 200f; // [g-C/L] (Bratbak 1985)

	// EPS
	protected static float specificMassEPS_WT; // [g-C/L]

	protected static float specificMassEPS_EpsMinus;

	// fraction of carbon that goes to EPS
	protected static float fEPS_WT;

	protected static float fEPS_EpsMinus;

	// Yield coefficients
	private static float Y = 0.44f; // [gX/gS]

	private static float YO = 2.66f; // [gO/gX]

	// Processes
	// Growth (biomass production)
	// protected static float uMaxWT = 0.39f; // [1/h] Trulear
	protected static float qSMax = 1.02f; // [1/h] from umax from Robinson

	// using yield from Beyenal

	// protected static float uMaxWT = 0.29f; // [1/h] Beyenal

	private static float KO = 1.18e-3f; // [g/L]

	// Computation parameters
	protected static float systemSize = 1000; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	// protected static float relativeMaximumRadius = 0.001f;
	protected static float relativeMaximumRadius = 4f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.5f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0; // NO DETACHMENT PRESENT IN THIS CASE

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumberWT = 50;

	protected static int initialParticleNumberEpsMinus = initialParticleNumberWT;

	protected static float finitshITerationTime = 480f;

	// program will only write to disk once, one hour before finishing
	// TODO change back
	protected static float timeForWritingToDisk = finitshITerationTime - 1;

	// protected static float timeForWritingToDisk = 1;
	protected static boolean runWithGraphics = false;

	// protected static float _mutationTime = finitshITerationTime*0.5f;
	protected static float mutationTime = 50f;

	protected static boolean mutationFlag = false;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// //1. Create the solutes
		// substrate
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// //2. Create the particulate species (solids)
		// X active mass EPS+
		ParticulateSpecies activeWT = new ParticulateSpecies("activeWT",
				specificMassX, Color.red);
		ParticulateSpecies activeEpsMinusMutant = new ParticulateSpecies(
				"activeEpsMinusMutant", specificMassX, Color.black);
		// a mutant strain that still produces EPS (for control)
		ParticulateSpecies activeControlMutant = new ParticulateSpecies(
				"activeControlMutant", specificMassX, Color.magenta);
		// EPS
		ParticulateSpecies eps_WT = new ParticulateSpecies("EPS_WT",
				specificMassEPS_WT, Color.yellow);
		// X active mass EPS-
		ParticulateSpecies activeEpsMinus = new ParticulateSpecies(
				"activeEpsMinus", specificMassX, Color.blue);
		// EPS
		ParticulateSpecies eps_EpsMinus = new ParticulateSpecies("EPS_EM",
				specificMassEPS_EpsMinus, Color.yellow);
		// //3. Create the biomass (bacterial) species
		// EPS -
		//
		ParticulateSpecies[] c1 = { activeEpsMinus, eps_EpsMinus };
		float fXVol_EPSMinus = 1.0f / (1.0f + fEPS_EpsMinus * specificMassX
				/ specificMassEPS_EpsMinus);
		float[] f1 = { 1.0f, 0 }; // use this to
		// start
		// particles
		BiomassSpecies epsMinus = new BiomassSpecies("epsMinus", c1, f1);
		epsMinus.setActiveMass(activeEpsMinus);
		epsMinus.setEpsMass(eps_EpsMinus);

		// EPS -
		ParticulateSpecies[] c2 = { activeWT, eps_WT, activeEpsMinusMutant,
				activeControlMutant };
		float fXVol_WT = 1.0f / (1.0f + fEPS_WT * specificMassX
				/ specificMassEPS_WT);
		float[] f2 = { 1.0f, 0, 0, 0 };
		BiomassSpecies wildType;
		if (mutationFlag) {
			wildType = new BiomassSpeciesForKeepingDistances("wildType", c2,
					f2, activeWT, activeEpsMinusMutant, mutationTime);
		} else {
			wildType = new BiomassSpeciesForKeepingDistances("wildType", c2,
					f2, activeWT, activeControlMutant, mutationTime);
		}
		wildType.setActiveMass(activeWT);
		wildType.setEpsMass(eps_WT);
		// //4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//
		// //5. Create the reactions
		// WT
		Reaction sUptakeWT = new Reaction("sUptakeWT", activeWT, qSMax, 1);
		sUptakeWT.addFactor(mO);
		// WT control mutant
		Reaction sUptakeControlMutant = new Reaction("sUptakeControlMutant",
				activeControlMutant, qSMax, 1);
		sUptakeControlMutant.addFactor(mO);
		// EPS -
		Reaction sUptakeEpsMinus = new Reaction("sUptakeEpsMinus",
				activeEpsMinus, qSMax, 1);
		sUptakeEpsMinus.addFactor(mO);
		// EPS - mutant
		Reaction sUptakeEpsMinusMutant = new Reaction("sUptakeEpsMinusMutant",
				activeEpsMinusMutant, qSMax, 1);
		sUptakeEpsMinusMutant.addFactor(mO);
		// // 6. Assign reaction to the species through ReactionStoichiometries
		// active mass EPS +
		NetReaction rsXactiveWT = new NetReaction(1);
		rsXactiveWT.addReaction(sUptakeWT, Y * 1 / (1 + fEPS_WT));
		activeWT.setProcesses(rsXactiveWT);
		// active mass EPS + control mutant
		NetReaction rsXactiveControlMutants = new NetReaction(1);
		rsXactiveControlMutants.addReaction(sUptakeControlMutant, Y * 1
				/ (1 + fEPS_WT));
		activeControlMutant.setProcesses(rsXactiveControlMutants);
		// EPS
		NetReaction rsEps_WT = new NetReaction(2);
		rsEps_WT.addReaction(sUptakeWT, Y * fEPS_WT / (1 + fEPS_WT));
		rsEps_WT.addReaction(sUptakeControlMutant, Y * fEPS_WT / (1 + fEPS_WT));
		eps_WT.setProcesses(rsEps_WT);
		// active mass EPS -
		NetReaction rsXactiveEpsMinus = new NetReaction(1);
		rsXactiveEpsMinus.addReaction(sUptakeEpsMinus, Y * 1
				/ (1 + fEPS_EpsMinus));
		activeEpsMinus.setProcesses(rsXactiveEpsMinus);
		// active mass EPS - mutant
		NetReaction rsXactiveEpsMinusMutant = new NetReaction(1);
		rsXactiveEpsMinusMutant.addReaction(sUptakeEpsMinusMutant, Y * 1
				/ (1 + fEPS_EpsMinus));
		activeEpsMinusMutant.setProcesses(rsXactiveEpsMinusMutant);
		// EPS
		NetReaction rsEps_EpsMinus = new NetReaction(1);
		rsEps_EpsMinus.addReaction(sUptakeEpsMinus, Y * fEPS_EpsMinus
				/ (1 + fEPS_EpsMinus));
		eps_EpsMinus.setProcesses(rsEps_EpsMinus);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(4);
		rsOxygen.addReaction(sUptakeWT, -YO * (1 - Y));
		rsOxygen.addReaction(sUptakeControlMutant, -YO * (1 - Y));
		rsOxygen.addReaction(sUptakeEpsMinus, -YO * (1 - Y));
		rsOxygen.addReaction(sUptakeEpsMinusMutant, -YO * (1 - Y));
		oxygen.setProcesses(rsOxygen);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(wildType);
		addBiomassSpecies(epsMinus);
		addSoluteSpecies(oxygen);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumberWT, initialParticleNumberEpsMinus };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * @param args
	 */
	public static void run() {
		// start
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new EpsProducersMutationKeepingDistances();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
		}
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			//
			Model.model().setCompulsoryTimeStep(timeForWritingToDisk);
			Model.model().setFinishIterationTime(finitshITerationTime);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			// writing only once
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new PovRayWriter()));
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			app.addStateWriter(new DistancesWriter());
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
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
			System.exit(-1);
		}
		try {
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}

	/**
	 * Set the variables and run using default output directory
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		runWithGraphics = true;
		// parse model parameters
		// fEPS_WT = 0.55f;
		// specificMassEPS_WT = specificMassX / 6f;
		// fEPS_EpsMinus = 1e-6f;
		// specificMassEPS_EpsMinus = specificMassX / 6f;
		//
		fEPS_WT = 0.5f;
		specificMassEPS_WT = specificMassX / 6f;
		// these are set for a EPS- mutant
		fEPS_EpsMinus = 1e-6f;
		specificMassEPS_EpsMinus = specificMassEPS_WT;
		// parse the initial frequencies
		int total = 100;
		initialParticleNumberWT = Math.round(total * 0.50f);
		initialParticleNumberEpsMinus = total - initialParticleNumberWT;
		//
		timeForWritingToDisk = 5f;
		Model.model().setSeed(100);
		//
		mutationTime = 50f;
		//
		run();
	}
}