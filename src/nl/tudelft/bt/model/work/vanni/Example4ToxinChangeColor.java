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
import nl.tudelft.bt.model.work.relatedness.Linear1;

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
public class Example4ToxinChangeColor extends ModelHandler {
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
	// protected static float substrateBulkConcentration = 1e-4f; // [g/L]

	// protected static float oxygenBulkConcentration = 8e-3f; // [g/L]

	// increasing oxygen to 80 not to make it limiting
	protected static float oxygenBulkConcentration = 8e1f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float oxygenDiffusivity = 8.33e6f; // [um^2/h]

	// Toxin (T) - Toxin bulk concentration (0 initial)-V.B.
	protected static float ToxinBulkConcentration = 0.0f; // [g/L]
	// Toxin diffusivity-V.B. Starting with same diffusivity as oxygen
	private static float ToxinDiffusivity = 0.03f ; //8.33e6f; // [um^2/h]
	
	private static float Toxinthreshold;

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 2.0e-7f; // [g/um^3]

	// Yield coefficients
	private static float Y = 0.44f; // [gX/gS]
	private static float YO = 2.66f; // [gO/gX]
	//-1/Y paper = Y/(Y0(1-Y))
	
	// Processes
	// Growth (biomass production)
	protected static float qsMax = 2.00f; // [1/h]   MUMAX = qsmax*Y = 0.88
	// protected static float uMax = 0.0547f; //[1/h]

	// Toxin production
	//protected static float qT = 0.1f; // [1/h]
	protected static float qT; //= 0.4f; // [1/h]
	// Toxin killing rate per Toxin unit concentration
	protected static float kT = 5;// [L/gT/hr]
	// Toxin loss rate per unit of susceptible species concentration
	protected static float BetaT = 0.0f; // [L/ gX /hr]
	protected static float alpha = 1.0f; // (gT/gX) Toxin to Biomass conversion
											// factor

	private static float KO = 1.18e-3f; // [g/L]

	// Computation parameters
	protected static float systemSize = 500/2f; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.01f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.5f;

	// other model parameters
	protected static int gridSide = 33; // mult0igrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0; // NO DETACHMENT PRESENT IN THIS CASE

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 200;

	// max cell mass
	protected static float cellmassmax;

	// min cell mass
	protected static float cellmassmin;

	// biofilm carachteristic length
	protected static float Lbac; 

	private static BiomassSpecies ToxinMinus;
	private static BiomassSpecies ToxinPlus;
	
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
		// Toxin species
		SoluteSpecies Toxin = new SoluteSpecies("Toxin", ToxinDiffusivity);
		// set up the simplest type of bulk concentration: constant
		Toxin.setBulkConcentration(new ConstantBulkConcentration(
				ToxinBulkConcentration));

		// 2. Create the particulate species (soliids)
		// X active mass
		ParticulateSpecies activeX = new ParticulateSpecies("activeX",
				specificMassX, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX = { activeX };
		float[] fractionalVolumeCompositionH1 = { 1.0f };

		// 2.1 V.B. Create the particulate species (soliids)
		// X2 active mass
		ParticulateSpecies activeX2 = new ParticulateSpecies("activeX2",
				specificMassX, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX2 = { activeX2 };
		float[] fractionalVolumeCompositionH2 = { 1.0f };

		// 3. Create the biomass species
		//ToxinMinus = new BiomassSpecies("ToxinMinus", spX,
		
		//		fractionalVolumeCompositionH1);
		
		//-----change color from Blue to Green if Toxin concentration is the above threshold----
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO);
		
		ToxinMinus = new BacteriocinBiomassSpecies("ToxinMinus", spX,
				fractionalVolumeCompositionH1, Toxin,
				Toxinthreshold, mO, Color.blue, Color.green);
		ToxinMinus.setActiveMass(activeX);

		// 3.1 V.B. Create the second biomass specie
		ToxinPlus = new BiomassSpecies("ToxinPlus", spX2,
				fractionalVolumeCompositionH2);
		ToxinPlus.setActiveMass(activeX2);

		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", activeX, qsMax, 1);
		growth.addFactor(mO);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		// similarly for activeX2
		Reaction growthX2 = new Reaction("growthX2", activeX2, qsMax, 1);
		growthX2.addFactor(mO);

		// Toxin production
		Reaction ToxinProduction = new Reaction("ToxinProduction", activeX2,
				qT, 1);

		// Toxin Killing effect
		Reaction ToxinKilling = new Reaction("ToxinKilling", activeX, kT, 1);
		ProcessFactor T1 = new Linear1(Toxin);
		ToxinKilling.addFactor(T1);

		// Toxin Depletion due to killing
		Reaction ToxinDepletion = new Reaction("ToxinDepletion", activeX,
				BetaT, 1);
		ProcessFactor T2 = new Linear1(Toxin);
		ToxinDepletion.addFactor(T2);

		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(2);
		rsXactive.addReaction(growth, Y);
		rsXactive.addReaction(ToxinKilling, -1);
		activeX.setProcesses(rsXactive);

		// similarly for activeX2 growth and investement in Toxin production
		NetReaction rsXactiveX2 = new NetReaction(2);
		rsXactiveX2.addReaction(growthX2, Y);
		rsXactiveX2.addReaction(ToxinProduction, -1);
		activeX2.setProcesses(rsXactiveX2);

		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(growth, -(YO * (1 - Y)));
		rsOxygen.addReaction(growthX2, -(YO * (1 - Y)));
		oxygen.setProcesses(rsOxygen);

		// Toxin
		NetReaction rsToxin = new NetReaction(2);
		rsToxin.addReaction(ToxinProduction, alpha);
		rsToxin.addReaction(ToxinDepletion, -1);
		Toxin.setProcesses(rsToxin);

		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(ToxinMinus);
		addBiomassSpecies(ToxinPlus);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(Toxin);
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
		int[] nCells = { initialParticleNumber, initialParticleNumber };
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

	// @ Override
	// public String toString() {
	// @Override public String toString() {
	// String s = new String();
	// s.concat("Lbac:"+ Lbac);
	// s.concat("alpha: " + alpha);
	// s.concat("diffusivity: " + ToxinDiffusivity);
	// s.concat("average MAX cell sixe : " + cellmassmax);
	// s.concat("kT : " + kT);
	// return s;
	// }

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		// read output directory from the command line
		if (args.length < 7) {
			throw new RuntimeException(
					"input arguments missing: \n"
							+ "1: output directory (CAUTION!!! directory will be erased \n"
							+ "2: seed for random number generator \n"
							+ "3: flag for running with graphics (1 on, 0 off) \n"
							+ "4: Kt \n"
							+ "5: Toxin Diffusivity \n"
							+ "6: 2D/3D \n"
							+ "7: qT \n");  
		}
		// parse inputs
		outputDirectory = args[0];
		parameterFile = outputDirectory +"Parameters.txt";
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		kT = Float.parseFloat(args[3]);
		ToxinDiffusivity=Float.parseFloat(args[4]);
		geometry = Integer.parseInt(args[5]);
		qT = Float.parseFloat(args[6]);
		
		cellmassmax = 3.1415f * (relativeMaximumRadius * systemSize)
		* (relativeMaximumRadius * systemSize) * systemSize / gridSide * specificMassX;
		
		cellmassmin = cellmassmax / 2.0f;
        
		Lbac = (alpha * kT * (cellmassmax + cellmassmin) / 2.0f)
		/ (2.0f * 3.1415f * ToxinDiffusivity);
		
		// Compute Toxin Threshold
		
		Toxinthreshold = qT/kT;
		
		//
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Example4ToxinChangeColor();
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
		System.out.println("Lbac is:" + Lbac);
		// System.out.println(Lbac);

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
			WriteParsToFile Parameters = new WriteParsToFile(parameterFile,
					true);
			// Parameters.writeToFile(app.toString());
			Parameters.writeToFile("Lbac: " + Lbac);
			Parameters.writeToFile("Alpha: " + alpha);
			Parameters.writeToFile("kT: " + kT);
			Parameters.writeToFile("Diffusivity: " + ToxinDiffusivity);
			Parameters.writeToFile("qT: " + qT);
			// /
			app.addTimedStateWriter(new PovRayWriter());
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
			spw.addSeries(new InverseOfDistancesSeries(ToxinMinus, ToxinPlus));
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
			 Model.model().setMaximumTimeStep(1);
			// Model.model().setFinishIterationTime(600);
			Model.model().setMaxHeight(150/1);
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