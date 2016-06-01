
package nl.tudelft.bt.model.work.carlos;

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
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
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
 * concentration of substrate.x`
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
public class Example1K extends TumorModelHandler {
	// All the model parameters are defined here as static attributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/Carlos/Documents/Results/Model/example1extend/";

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
	//protected static float substrateBulkConcentration = 1e-4f; // [g/L]
	//glucose plasma concentration = 5.5mM = 0.99g/L i.e. 1g/L
	protected static float substrateBulkConcentration = 1f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	// Glucose diffusivity (Casciari et al 1988) 1.1e-6 cm2/s = 3.96e5um2/h mouse; 5.5e-7-2.3e-7 cm2/s = 1.98e4um2/h to 8.2e3 um2/h
	private static float substrateDiffusivity = 8.2e3f ; // [um^2/h]

	// 

	//
	// Particulate species (biomass X)
	protected static float specificMassX1 = 100f; // [g/L]
	protected static float specificMassX2 = 100f; // [g/L]


	// Yield coefficients
	private static float YXS = 0.045f; // [gX/gS]

	// Processes
	// Growth (biomass production)
	protected static float uMax1 = 0.01f; // [1/h]
	protected static float uMax2 = 0.001f; //[1/h]

	private static float KS = 3.5e-4f; // [g/L]
	
	// Death rates
	
	protected static float kd1 = 0.001f;
	protected static float kd2 = 0.01f;

	// Computation parameters 
	protected static float systemSize = 2000; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 5f/systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.05f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 1f/systemSize;

	// other model parameters
	protected static int gridSide = 129; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side
	
	protected static int numberofCellGroups = 2;


	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 1e-5f; // 

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 500;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		// set up the simplest type of bulk concentration: constant
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));

		// 2. Create the particulate species (solids)
		// X active mass
		ParticulateSpecies activeX1 = new ParticulateSpecies("activeX1",
				specificMassX1, Color.green);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX1 = { activeX1 };
		float[] fractionalVolumeCompositionH1 = { 1.0f };
		
		ParticulateSpecies activeX2 = new ParticulateSpecies("activeX2",
				specificMassX2, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX2 = { activeX2 };
		float[] fractionalVolumeCompositionH2 = { 1.0f };

		// 3. Create the biomass species
		BiomassSpecies speciesX1 = new BiomassSpecies("speciesX1", spX1,
				fractionalVolumeCompositionH1);
		speciesX1.setActiveMass(activeX1);
		//speciesX.getColorFromGrowth();

		BiomassSpecies speciesX2 = new BiomassSpecies("speciesX2", spX2,
				fractionalVolumeCompositionH2);
		speciesX2.setActiveMass(activeX2);
		//speciesY.getColorFromGrowth();
	
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growthX1 = new Reaction("growthX1", activeX1, uMax1, 1);
		growthX1.addFactor(mS);
		
		Reaction growthX2 = new Reaction("growthX2", activeX2, uMax2, 1);
		growthX2.addFactor(mS);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		
		
		// Endogenous Decay
		Reaction BiomassDecayX1 = new Reaction("BiomassDecayX1", activeX1, kd1, 0);
		Reaction BiomassDecayX2 = new Reaction("BiomassDecayX2", activeX2, kd2, 0);


		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactiveX1 = new NetReaction(2);
		rsXactiveX1.addReaction(growthX1, 1);
		rsXactiveX1.addReaction(BiomassDecayX1, -1);
		activeX1.setProcesses(rsXactiveX1);
		
		NetReaction rsXactiveX2 = new NetReaction(2);
		rsXactiveX2.addReaction(growthX2, 1);
		rsXactiveX2.addReaction(BiomassDecayX2, -1);
		activeX2.setProcesses(rsXactiveX2);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsGlucose = new NetReaction(2);
		rsGlucose.addReaction(growthX1, -(1 / YXS));
		rsGlucose.addReaction(growthX2, -(1 / YXS));
		substrate.setProcesses(rsGlucose);
		
	
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesX1);
		addBiomassSpecies(speciesX2);
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
		int[] nCells = { initialParticleNumber/2,initialParticleNumber*2 };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc) 
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(rdetach);
		setDetachmentHandler(df);
		Model.model().getBiomassContainer().turnSloughingOff();
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Example1K();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries runLengthX = new RunLengthXSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		app = new SeriesVizualizer(app, runLengthX);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
            app.addTimedStateWriter(new ImageWriter());
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(runLengthX);
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