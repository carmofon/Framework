package nl.tudelft.bt.model.examplesroland;

import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;

/**
 * Theis program runs the model of competition between Heterotrophic and
 * autotrophic organisms as described in Morgenroth and Wilderer (2000)
 * 
 * @author Joao and Roland
 */
public class AutotrophsHeterotrophs extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	//output directory name (
	protected static String outputDirectory = "/Users/jxavier/results/AutotrophsHeterotrophs";
	//protected static String outputDirectory = "D:\\0_Delft\\results\\test2";

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

	//Solute species
	//Substrate (S) - the only solute species used here
	protected static float oxygenBulkConcentration = 8f; //[gO/m3]

	protected static float substrateBulkConcentration = 10f; //[gCOD-S/m3]

	protected static float ammoniumBulkConcentration = 13f; //[gN-NH4/m3]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float oxygenDiffusivity = 2.1e-4f; //[m^2/d]

	private static float substrateDiffusivity = 5.8e-5f; //[m^2/d]

	// 

	//
	//Particulate species (biomass X)
	protected static float specificMassH = 70000f; //[gCOD-H/m3]

	protected static float specificMassA = 70000f; //[gCOD-A/m3]
	protected static float specificMassEPS = 10000f; //[gCOD-EPS/m3]

	//Yield coefficients
	private static float YH = 0.4f; //[gCOD-H/gCOD-S]

	private static float YA = 0.22f; //[gCOD-A/gCOD-N]
	private static float YEPS = 0.3f; //[gCOD-EPS/gCOD-H]

	// Processes
	//Growth (biomass production)
	private static float uMaxH = 4.8f; //[1/day]

	private static float uMaxA = 0.95f; //[1/day]

	private static float kH = 0.1f; //[1/day]

	private static float kA = 0.1f; //[1/day]

	private static float bH = 0.2f; //[1/day]
	private static float bA = 0.05f; //[1/day]

	private static float KO2 = 0.35f; //[gO/m3]

	private static float KS = 5.0f; //[gCOD/m3]

	private static float KNH4 = 1.0f; //[gN-NH4/m3]

	// Computation parameters
	protected static float systemSize = 2e-3f; // [m]

	//relativeMaximumRadius defines the maximum radius of the biomass particles
	//in relation to the system size
	//the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.006f;

	//Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 1e7f; // 

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 10;

	//reactor properties
	private static float residenceTime = 1f; // [day]
	private static float reactorVolume = 1f; // [m3]
	private static float carrierArea = 200f / reactorVolume; // [m2]

	// the compulsory time step
	private static float outputEvery = 0.5f;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		//1. Create the solutes
		//OXYGEN
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		//set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		//SUBSTRATE
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		//set up the simplest type of bulk concentration: constant
		//		substrate.setBulkConcentration(new ConstantBulkConcentration(
		//				substrateBulkConcentration));
		substrate.setBulkConcentration(new DynamicBulkConcentrationImplicit(
				substrateBulkConcentration));
		//AMMONIUM
		SoluteSpecies ammonium = new SoluteSpecies("ammonium",
				substrateDiffusivity);
		//set up the simplest type of bulk concentration: constant
		ammonium.setBulkConcentration(new ConstantBulkConcentration(
				ammoniumBulkConcentration));
		//2. Create the particulate species (solids)
		//H active mass of heterotrophs
		ParticulateSpecies activeH = new ParticulateSpecies("activeH",
				specificMassH, Color.red);
		ParticulateSpecies inertH = new ParticulateSpecies("inertH",
				specificMassH, Color.red);
		ParticulateSpecies epsH = new ParticulateSpecies("EPS",
				specificMassEPS, Color.green);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spH = {activeH, inertH, epsH};
		float[] fractionalVolumeCompositionH1 = {1.0f, 0f, 0f};
		//3. Create the biomass species
		BiomassSpecies speciesH = new BiomassSpecies("speciesH", spH,
				fractionalVolumeCompositionH1);
		//the following are set only for colour purposes
		speciesH.setActiveMass(activeH);
		speciesH.setInertMass(inertH);
		speciesH.setEpsMass(epsH);
		//A active mass of autotrophs
		ParticulateSpecies activeA = new ParticulateSpecies("activeA",
				specificMassH, Color.blue);
		ParticulateSpecies inertA = new ParticulateSpecies("inertA",
				specificMassH, Color.blue);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spA = {activeA, inertA};
		float[] fractionalVolumeCompositionA = {1.0f, 0f};
		//3. Create the biomass species
		BiomassSpecies speciesA = new BiomassSpecies("speciesA", spA,
				fractionalVolumeCompositionA);
		speciesA.setActiveMass(activeA);
		speciesA.setInertMass(inertA);
		//4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO2);
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor mNH4 = new Saturation(ammonium, KNH4);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//
		//5. Create the reactions
		//growth of heterotrophs
		Reaction growthH = new Reaction("growthH", activeH, uMaxH, 2);
		growthH.addFactor(mO);
		growthH.addFactor(mS);
		//inactivation of heterotrophs
		Reaction inactivationH = new Reaction("inactivationH", activeH, kH, 0);
		//endogenous respiration of heterotrophs
		Reaction respirationH = new Reaction("respirationH", activeH, bH, 1);
		respirationH.addFactor(mO);
		//growth of autotrophs
		Reaction growthA = new Reaction("growthA", activeA, uMaxA, 2);
		growthA.addFactor(mO);
		growthA.addFactor(mNH4);
		//inactivation of autotrophs
		Reaction inactivationA = new Reaction("inactivationA", activeA, kA, 0);
		//endogenous respiration of autotrophs
		Reaction respirationA = new Reaction("respirationA", activeA, bA, 1);
		respirationA.addFactor(mO);
		//
		//6. Assign reaction to the species through ReactionStoichiometries
		//active mass of H
		NetReaction rsActiveH = new NetReaction(3);
		rsActiveH.addReaction(growthH, 1);
		rsActiveH.addReaction(inactivationH, -1);
		rsActiveH.addReaction(respirationH, -1);
		activeH.setProcesses(rsActiveH);
		//inert mass of H
		NetReaction rsInertH = new NetReaction(1);
		rsInertH.addReaction(inactivationH, 1);
		inertH.setProcesses(rsInertH);
		//EPS of H
		NetReaction rsEpsH = new NetReaction(1);
		rsEpsH.addReaction(growthH, YEPS);
		epsH.setProcesses(rsEpsH);
		//active mass of A
		NetReaction rsActiveA = new NetReaction(3);
		rsActiveA.addReaction(growthA, 1);
		rsActiveA.addReaction(inactivationA, -1);
		rsActiveA.addReaction(respirationA, -1);
		activeA.setProcesses(rsActiveA);
		//inert mass of A
		NetReaction rsInertA = new NetReaction(1);
		rsInertA.addReaction(inactivationA, 1);
		inertA.setProcesses(rsInertA);
		// This defines that biomass growth rate is 1*rX
		//
		//assign reaction stoichiometry to the solutes
		//substrate
		NetReaction rsOxygen = new NetReaction(4);
		rsOxygen.addReaction(growthH, -1 / YH + 1 + YEPS);
		rsOxygen.addReaction(growthA, -(4.57f - YA) / YA);
		rsOxygen.addReaction(respirationH, -1);
		rsOxygen.addReaction(respirationA, -1);
		oxygen.setProcesses(rsOxygen);
		//SUBSTRATE
		NetReaction rsSubstrate = new NetReaction(1);
		rsSubstrate.addReaction(growthH, -1 / YH);
		substrate.setProcesses(rsSubstrate);
		//AMMONIUM
		NetReaction rsAmmonium = new NetReaction(1);
		rsAmmonium.addReaction(growthA, -1 / YA);
		ammonium.setProcesses(rsAmmonium);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		//7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesH);
		addBiomassSpecies(speciesA);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(substrate);
		addSoluteSpecies(ammonium);
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
		int[] nCells = {initialParticleNumber, initialParticleNumber};
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
		MultigridVariable.setSteps(100, 1000);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new AutotrophsHeterotrophs();
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
		// add vizualizer for solutes rates
		app = new BulkConcentrationVizualizer(app);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			//initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new PovRayWriter()));
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addTimedStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			//app.addStateWritter(new DetachmentLevelSetWriter());
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
			// set the reactor dimensions
			app.setReactorParameters(residenceTime, carrierArea, reactorVolume);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// set the times to write out
			Model.model().setCompulsoryTimeStep(outputEvery);
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}