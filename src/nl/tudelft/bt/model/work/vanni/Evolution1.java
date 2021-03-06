package nl.tudelft.bt.model.work.vanni;

import java.awt.Color;
import java.io.File;
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
import nl.tudelft.bt.model.util.ColorMaps;
import nl.tudelft.bt.model.work.relatedness.Linear1;
import nl.tudelft.bt.model.work.vanni.MutatorBiomassSpecies;

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
public class Evolution1 extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/bucciv/results/test8/";
	protected static String parameterFile = outputDirectory +"Parameters.txt";

	// WARNING: the contents of the output directory will be deleted!!
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
	protected static float oxygenBulkConcentration ; //= 8e1f; // [gO/um^3]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h 8.33e6
	private static float oxygenDiffusivity ;//= 4e4f; // [um^2/h]

	// Particulate species (biomass X)
	protected static float specificMassX = 1.5e-7f; // [g/um^3]
	//protected static float specificMassX = 2.0e-7f; // [g/um^3]
	
	// Yield coefficients
	private static float Y = 0.5f; // [gX/gX]
	private static float YO = 0.5f; // [gO/gX]

	// Processes
	// Growth (biomass production)
	protected static float qsMax = 2f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]

	private static float KO = 3.5e-14f; // [gO/um^3]
	// Computation parameters
	protected static float systemSize = 250f; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 1.0f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um 
	protected static float relativeBoundaryLayer =  25f / systemSize;

	// other model parameters
	protected static int gridSide = 65; // mult0igrid grid side
	// protected static int gridSide = 33; // multigrid grid side
	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	// Detachment h_star is the steady state biofilm height 
	protected static float h_star; //[um]
	
	protected static float rdetach; // = Y*qsMax*specificMassX / h_star; // 0 = NO DETACHMENT PRESENT IN THIS CASE
 
	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 130;
	
	protected static float InitialRatio;

	// max cell mass
	protected static float cellmassmax;

	// min cell mass
	protected static float cellmassmin;

	protected static double Delta;
	protected static float mutationRate;
	
	private static EvolutionaryBiomassSpecies BiomassX;
	//private static BiomassSpecies ToxinMinus;
	//private static BiomassSpecies ToxinPlus;
	
	private static MutatorBiomassSpecies speciesX;
	
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// Oxygen species
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		
		// 2. Create the particulate species (soliids)
		// X active mass
		ParticulateSpecies activeX = new ParticulateSpecies("activeX",
				specificMassX, Color.blue);
		
		//ParticulateSpecies P1 = new ParticulateSpecies("P1", 1.0f, Color.red );
		
		//ParticulateSpecies[] spX = {activeX, P1};
		ParticulateSpecies[] spX = {activeX};
		//float[] fractionalVolumeCompositionH1 = { 0.9f, 0.1f };
		float[] fractionalVolumeCompositionH1 = { 1.0f};
		// 3. Create the biomass species
		BiomassX = new EvolutionaryBiomassSpecies("BiomassX", spX,
				fractionalVolumeCompositionH1, 0.1f);
		BiomassX.setActiveMass(activeX);
	    //BiomassX.setParameterValue(P1);	

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO);
		//Reaction P1Constant = new Reaction("P1Constant", P1, 0, 1);
		//ProcessFactor _P1 = new Linear1(P1);
		// The Saturation class creates; a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", activeX, qsMax, 2);
		growth.addFactor(mO);
		//growth.addFactor(_P1);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(1);
		rsXactive.addReaction(growth, Y);
		activeX.setProcesses(rsXactive);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(growth, -(YO * (1-Y)));
		oxygen.setProcesses(rsOxygen);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		// EPS stochiometry
		//NetReaction rsP1 = new NetReaction(1);
		//rsP1.addReaction(P1Constant, 1);
		//P1.setProcesses(rsP1);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(BiomassX);
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
		//int[] nCells = { (int) ((1-InitialRatio)*initialParticleNumber), (int) (InitialRatio*initialParticleNumber) };
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
	 * @throws IOException
	 */
	
	
	
	public static void main(String[] args) throws IOException {
		// read output directory from the command line
		if (args.length < 12) {
			throw new RuntimeException(
					"input arguments missing: \n"
							+ "1: output directory (CAUTION!!! directory will be erased \n"
							+ "2: seed for random number generator \n"
							+ "3: flag for running with graphics (1 on, 0 off) \n"
							+ "4: qSmax \n"
							+ "5: Y \n"
							+ "6: 2D/3D \n"
							+ "7: Y0 f \n"
							+ "8: Oxygen Bulk Concentration \n"
							+ "9: Oxygen Diffusivity \n"
							+ "10: K0"
							+ "11: Mutation Rate \n"
							+ "12: Steady State height h*"	);
		}
		// parse inputs
		outputDirectory = args[0];
		parameterFile = outputDirectory +"Parameters.txt";
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		qsMax = Float.parseFloat(args[3]);
		Y=Float.parseFloat(args[4]);
		geometry = Integer.parseInt(args[5]);
		YO = Float.parseFloat(args[6]);
		oxygenBulkConcentration = Float.parseFloat(args[7]);
		oxygenDiffusivity = Float.parseFloat(args[8]);
		KO =  Float.parseFloat(args[9]);
		mutationRate = Float.parseFloat(args[10]);
		h_star = Float.parseFloat(args[11]);
		
		rdetach = Y*qsMax*specificMassX / h_star; // 0 = NO DETACHMENT PRESENT IN THIS CASE
	
		cellmassmax = 3.1415f * (relativeMaximumRadius * systemSize)
		* (relativeMaximumRadius * systemSize) * systemSize / gridSide * specificMassX;
		
		cellmassmin = cellmassmax / 2.0f;
		
		Float delta_square = new Float (oxygenBulkConcentration * oxygenDiffusivity * Y/(YO*(1-Y))) / ((qsMax * Y) * specificMassX * (relativeBoundaryLayer * systemSize) * (relativeBoundaryLayer * systemSize));
		Delta = Math.sqrt(delta_square.doubleValue());
		//
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Evolution1();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
		}
		// print Lbac
		// System.out.println(Lbac);
		System.out.println("Delta is:" + Delta);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			// write initial parameters to file
			//WriteParsToFile Parameters = new WriteParsToFile(parameterFile,
			//		true);
			// Parameters.writeToFile(app.toString());
			//Parameters.writeToFile("Lbac: " + Lbac);
			//Parameters.writeToFile("Alpha: " + alpha);
			//Parameters.writeToFile("kT: " + kT);
			//Parameters.writeToFile("Diffusivity: " + ToxinDiffusivity);
			//Parameters.writeToFile("h_star: " + h_star);
			//Parameters.writeToFile("Delta: " + Delta);
			// /
			//app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
		    app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			//app.addStateWritter(new DetachmentLevelSetWriter());
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// the simulation parameters writer
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
			//spw.addSeries(new InverseOfDistancesSeries(ToxinMinus, ToxinPlus));
			app.addStateWriter(spw);
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
			//Model.model().setCompulsoryTimeStep(240f);
			//Model.model().setFinishIterationTime(600);
			//Model.model().setMaxThickness(150);
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