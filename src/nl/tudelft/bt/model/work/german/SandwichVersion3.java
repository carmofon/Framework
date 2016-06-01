package nl.tudelft.bt.model.work.german;

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
 * Example 1: Growth of a monospecies biofilm with a single substrate species.
 * No detachment forces present and growth is represented in 2D space. The
 * biofilm developed shows a rough morphology (heterogeneous shape) as a
 * consequance of high diffusion limitation. The effect of diffusion limitation
 * is represented by the dimesionless G number, defined as
 * 
 * G = Ly^2*u_max*C_X_max/(D_S*C_S_bulk) where: Ly - size of the system u_max -
 * maximum specific growth rate of microorganisms C_X_max - Specific mass of
 * biomass D_S - Diffusion coefficient of substrate C_S_bulk - bulk
 * concentration of substrate.
 * 
 * Parameters taken from
 * 
 * Picioreanu, C., Van Loosdrecht, M. C. M. and Heijnen, J. (1998) Mathematical
 * modelling of biofilm structure with a hybrid differential-discrete cellular
 * automaton approach. Biotech Bioeng, 58, 101-116.
 * 
 * Note: um representes micro-meter
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SandwichVersion3 extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/cbioguest/results/test1/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "e:/results/example1/" directory "e:\results" must
	// exist in your computer. The subdirectory "example1" will be created
	// if it is non-existent. If it exists, its contents will be deleted
	// during the program initialization

	// geometry (default is 2D) - change value to 3 for 3D
	protected static int geometry = 2;

	// Solute species
	// Substrate (S) - the only solute species used here
	//protected static float substrateBulkConcentration = 1e-4f; // [g/L]
	protected static float substrateBulkConcentration = 1e-4f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float substrateDiffusivity = 5.76e6f; // [um^2/h]
	// 
	protected static SoluteSpecies substrate;
	//
	// Particulate species (biomass X)
	protected static float specificMassX = 70f; // [g/L]

	// Yield coefficients
	private static float YXS = 0.045f; // [gX/gS]

	// Processes
	// Growth (biomass production)
	protected static float uMax = 0.1f; // [1/h]
	protected static float uMaxY = 0.1f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]

	private static float KS = 3.5e-4f; // [g/L]
	//private static float KSY = 3.5e-4f; // [g/L]

	private static float kDeathX = 0.001f; // [h-1]
	//private static float kDeathY = 0.001f; // [h-1]

	
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
	protected static int initialParticleNumber = 10;



	
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// substrate
		substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		// set up the simplest type of bulk concentration: constant
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));

		// 2. Create the particulate species (solids)
		// X active mass
		ParticulateSpecies activeX = new ParticulateSpecies("activeX",
				specificMassX, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX = { activeX };
		float[] fractionalVolumeCompositionH1 = { 1.0f };

		// 3. Create the biomass species
		BiomassSpecies speciesX = new BiomassSpecies("speciesX", spX,
				fractionalVolumeCompositionH1);
		speciesX.setActiveMass(activeX);
		//speciesX.getColorFromGrowth();

		// Create the second particulate species (solids)
		// Y active mass
		//ParticulateSpecies activeY = new ParticulateSpecies("activeY",
		//		specificMassX, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		//ParticulateSpecies[] spY = { activeY };
		//float[] fractionalVolumeCompositionY = { 1.0f };

		// Create the second biomass species
		//BiomassSpecies speciesY = new BiomassSpecies("speciesY", spY,
		//		fractionalVolumeCompositionY);
		//speciesY.setActiveMass(activeY);
		//speciesY.getColorFromGrowth();
		
		
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", activeX, uMax, 1);
		growth.addFactor(mS);

		Reaction deathX = new Reaction("deathX", activeX, kDeathX, 0);
		
		//Reaction growthY = new Reaction("growthY", activeY, uMaxY, 1);
		//growthY.addFactor(mS);

		//Reaction deathY = new Reaction("deathY", activeY, kDeathY, 0);

		
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(2);
		rsXactive.addReaction(growth, 1);
		rsXactive.addReaction(deathX, -1);
		activeX.setProcesses(rsXactive);

		
		//NetReaction rsYactive = new NetReaction(2);
		//rsYactive.addReaction(growthY, 1);
		//rsYactive.addReaction(deathY, -1);
		//activeY.setProcesses(rsYactive);

		
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(growth, -(1 / YXS));
		//rsSubstrate.addReaction(growthY, -(1 / YXS));
		substrate.setProcesses(rsSubstrate);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesX);
		//addBiomassSpecies(speciesY);
		addSoluteSpecies(substrate);
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
		MultigridVariable.setSteps(5, 50);
		// create a handle for the application, which will be decorated
		ApplicationComponent app = new SandwichVersion3();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no visuals are desired
		// start decorating the application
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
			app.addTimedStateWriter(new ImageWriter(substrate));
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