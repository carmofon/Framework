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
 * This program should predict the structural development of an autotrophic
 * biofilm. At first it includes only growth of bacteria on oxygen and ammonium,
 * then as much processes are added until the structural development is alike
 * the one described in Staudt et al. 2004 Kinetic parameters are taken from
 * Horn and Hempel 1997 Growth and decay in an auto-/heterotrophic biofilm, Wat.
 * Res. 31(9) 1997.
 * 
 * In this second version it should be tried do describe the cultivation "Amm7",
 * conducted by Christian Staudt. To figure out the parameters the data of this
 * cultivation have been used for a mass balance. From here a new yield could be
 * calculated for EPS+Bacteria-production on consumption of ammonia consumed.
 * The ration between bacteria and EPS is in the range between 1.5 and 20 as can
 * be derived from the distribution profiles of EPS and bacteria in "Amm9".
 * 
 * @author Roland, Software written by Joao Xavier
 */
public class AutotrophsEps2 extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	//output directory name (
	//protected static String outputDirectory =
	// "D:\\0_Delft\\results\\AutoEps";
	protected static String outputDirectory = "c:\\roland\\results\\AutoEps";

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

	//SOLUTES
	//Bulk concentrations
	protected static float oxygenBulkConcentration = 8f; //[gO/m3]

	//	the value for the bulk concentration is assumed to be 30g/m^3, actually
	// it should
	//	be 45.1g/m^3, but due to numeric problems and experimental data showing
	// until day
	//	60 a sufficient supply of ammonium, a constant concentration is assumed

	protected static float ammoniumBulkConcentration = 45.1f; //[gN-NH4/m3]

	// Diffusivities
	private static float oxygenDiffusivity = 2.1e-4f; //[m^2/d]

	private static float ammoniumDiffusivity = 1.8e-4f; //[m^2/d]

	//SOLIDS
	//Particulate species (biomass X)
	//in AQUASIM both: 250kg/m^3, but eps_s=0.06
	//now assumed: BM_ges=341kg/m^3*0.06
	protected static float specificMassA = 24800f; //[g COD-A/m3]

	protected static float specificMassEpsA = 24800f; //[gCOD-EPS/m3]

	//Yield coefficients
	//estimated from data of C. Staudt (YEpsA from Amm7, YEps from Amm9)
	//??Is the value YEpsA for NH4 or N-NH4??
	private static float YEpsA = 0.77f; //[gCOD-A/gCOD-N-NH4]

	private static float YEps = 1.5f; //[gCOD-EPS/gCOD-A]

	//PROCESSES
	//Growth (biomass production)
	private static float uMaxA = 0.14f; //[1/day]

	private static float KO2 = 0.5f; //[gO/m3]

	private static float KNH4 = 0.5f; //[gN-NH4/m3]

	//COMPUTATIONAL PARAMETERS
	protected static float systemSize = 2e-4f; // [m]

	//relativeMaximumRadius defines the maximum radius of the biomass particles
	//in relation to the system size
	//the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 1.2e-6f / systemSize;

	//Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0f / systemSize;

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0f; // 

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 10;

	//reactor properties
	private static float residenceTime = 0.25f; // [day]

	private static float reactorVolume = 1.455e-3f; // [m3]

	private static float carrierArea = 0.225f; // [m2]

	//cut off biofilm if following height is reached:
	private static float BiofilmCutOffHeight = 1.028e-4f;

	// the compulsory time step
	private static float outputEvery = 0.5f;

	//at which time shall the simulation stop?
	private static float endSimulationTime = 60f;

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

		//AMMONIUM
		SoluteSpecies ammonium = new SoluteSpecies("ammonium",
				ammoniumDiffusivity);
		//set up the simplest type of bulk concentration: constant
		ammonium.setBulkConcentration(new ConstantBulkConcentration(
				ammoniumBulkConcentration));
		//		ammonium.setBulkConcentration(new DynamicBulkConcentrationImplicit(
		//		ammoniumBulkConcentration));

		//2. Create the particulate species (solids)
		//A active mass of autotrophs
		ParticulateSpecies activeA = new ParticulateSpecies("activeA",
				specificMassA, Color.red);
		ParticulateSpecies EpsA = new ParticulateSpecies("EPS",
				specificMassEpsA, Color.green);
		//		ParticulateSpecies inertA = new ParticulateSpecies("inertA",
		//				specificMassH, Color.blue);

		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spA = { activeA, EpsA };
		float[] fractionalVolumeCompositionA = { 1.0f, 0f };
		//3. Create the biomass species
		BiomassSpecies speciesA = new BiomassSpecies("speciesA", spA,
				fractionalVolumeCompositionA);
		speciesA.setActiveMass(activeA);
		speciesA.setEpsMass(EpsA);
		//
		//4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO2);
		ProcessFactor mNH4 = new Saturation(ammonium, KNH4);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//
		//5. Create the reactions
		//growth of autotrophs
		Reaction growthA = new Reaction("growthA", activeA, uMaxA, 2);
		growthA.addFactor(mO);
		growthA.addFactor(mNH4);
		//
		//6. Assign reaction to the species through ReactionStoichiometries
		//active mass of A
		NetReaction rsActiveA = new NetReaction(1);
		rsActiveA.addReaction(growthA, 1f);
		activeA.setProcesses(rsActiveA);
		//EPS of A
		NetReaction rsEpsA = new NetReaction(1);
		rsEpsA.addReaction(growthA, YEps);
		EpsA.setProcesses(rsEpsA);

		//assign reaction stoichiometry to the solutes
		//substrate
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(growthA, 1f + YEps - 4.4f * ((1f + YEps) / YEpsA));
		oxygen.setProcesses(rsOxygen);
		//AMMONIUM
		NetReaction rsAmmonium = new NetReaction(1);
		rsAmmonium.addReaction(growthA, -(1f + YEpsA) / YEpsA);
		ammonium.setProcesses(rsAmmonium);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		//7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesA);
		addSoluteSpecies(oxygen);
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
		int[] nCells = { initialParticleNumber };
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
		System.out.println("Starting");
		System.out.println("args[0] = " + args[0]);
		// check input arguments
		boolean graficsOn = true;
		if (args.length > 0) {
			if (args[0].equals("-nodisplay")) {
				graficsOn = false;
				System.out.println("Running without graphics");
			}
			else
				throw new RuntimeException("Illegal argument " + args[0]);
		}
		//
		MultigridVariable.setSteps(100, 1000);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new AutotrophsEps2();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		if (graficsOn) {
			// The following code will be omitted if no vizuals are desired
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
			// add vizualizer for bulk concentration of solutes
			app = new BulkConcentrationVizualizer(app);
			// add vizualizer for solutes rates
			//		app = new SoluteRateSeriesVizualizer(app);
			//add buttons
			app = new VizualModelControler(app);
		}
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
			//Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(endSimulationTime);
			Model.model().setVerticalCutoffSize(BiofilmCutOffHeight);
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}