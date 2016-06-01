package nl.tudelft.bt.model.work.will;

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
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * Tumor with two neutral cell types
 * For now X will be producer and Y will be consumer
 * 
 * @author Will Chang
 */
public class TwoCellTypesDecaySignal extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/wkc/xavierlab/results/tests";

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
	// Oxygen (O)
	protected static float oxygenBulkConcentration = 1e-4f; // [g/L] produces smooth tumor
	//protected static float oxygenBulkConcentration = 1e-20f; // [g/L] produces fingers
	
	// Signal (S)
	protected static float signalBulkConcentration = 1e-6f;

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float oxygenDiffusivity = 5.76e6f; // [um^2/h]
	private static float signalDiffusivity = 1.76e6f; // [um^2/h]

	// 

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 3000f; // [g/L]

	// Yield coefficients
	// oxygen
	private static float YXO = 0.045f; // [gX/gS]
	// growth signal
	private static float YXS = 0.45f;

	// Processes
	// Growth (biomass production)
	protected static float uMax = 0.1f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]
	
	// Decay (biomass destruction)
	protected static float kd = 1e-5f;
	
	// Growth signal production rate
	protected static float Rs = 10f;
	
	// Growth signal production cost (this should be < 0)
	protected static float c = -0.01f;

	// Half-saturation constant for oxygen growth
	private static float KO = 3.5e-4f; // [g/L]
	
	// Half-saturation constant for growth signal growth
	private static float KS = 3.5e-4f; //[g/L]

	// Computation parameters
	protected static float systemSize = 100; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.005f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.05f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0; // NO DETACHMENT PRESENT IN THIS CASE

	// initial number of particles in the system (inoculum)
	protected static int initialCellNumber = 10;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen",
				oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// growth signal
		SoluteSpecies signal = new SoluteSpecies("growth signal", signalDiffusivity);
		signal.setBulkConcentration(new ConstantBulkConcentration(signalBulkConcentration));

		// 2a. Create the particulate species (solids)
		// X active mass
		ParticulateSpecies activeX = new ParticulateSpecies("activeX",
				specificMassX, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX = { activeX };
		float[] fractionalVolumeCompositionH1 = { 1.0f };


		// 2b. Create the particulate species (solids)
		// Y active mass
		ParticulateSpecies activeY = new ParticulateSpecies("activeY",
				specificMassX, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spY = { activeY };
		float[] fractionalVolumeCompositionY = { 1.0f };

		// 3a. Create the biomass speciesX
		BiomassSpecies speciesX = new BiomassSpecies("speciesX", spX,
				fractionalVolumeCompositionH1);
		speciesX.setActiveMass(activeX);

		// 3b. Create the biomass speciesY
		BiomassSpecies speciesY = new BiomassSpecies("speciesY", spY,
				fractionalVolumeCompositionY);
		speciesY.setActiveMass(activeY);
		
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(oxygen, KO);
		ProcessFactor signalS = new Saturation(signal, KS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growthX = new Reaction("growthX", activeX, uMax, 2);
		growthX.addFactor(mS);
		growthX.addFactor(signalS);

		Reaction growthY = new Reaction("growthY", activeY, uMax, 2);
		growthY.addFactor(mS);
		growthY.addFactor(signalS);
		
		// decay
		Reaction decayX = new Reaction("decayX", activeX, kd, 0);
		Reaction decayY = new Reaction("decayY", activeY, kd, 0);
		
		// growth signal production
		Reaction produceSignal = new Reaction("produceSignal", activeX, Rs * uMax, 2);
		produceSignal.addFactor(mS);
		produceSignal.addFactor(signalS);

		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(3);
		rsXactive.addReaction(growthX, 1);
		rsXactive.addReaction(decayX, -1);
		rsXactive.addReaction(produceSignal, c);
		activeX.setProcesses(rsXactive);
		//
		NetReaction rsYactive = new NetReaction(2);
		rsYactive.addReaction(growthY, 1);
		rsYactive.addReaction(decayY, -1);
		activeY.setProcesses(rsYactive);
		//
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(growthX, -(1 / YXO));
		rsOxygen.addReaction(growthY, -(1 / YXO));
		oxygen.setProcesses(rsOxygen);
		// This defines that oxygen consumption rate is -(1 / YXO)*(rX + rY)
		// growth signal
		NetReaction rsSignal = new NetReaction(3);
		rsSignal.addReaction(growthX, -(1 / YXS));
		rsSignal.addReaction(growthY, -(1 / YXS));
		rsSignal.addReaction(produceSignal, 1);
		signal.setProcesses(rsSignal);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesX);
		addBiomassSpecies(speciesY);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(signal);
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
		int[] nCells = { initialCellNumber, initialCellNumber };
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
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TwoCellTypesDecaySignal();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		app = new SeriesVizualizer(app, thickness);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
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
}
