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
public class Example9ToxinDiversity extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/bucciv/results/test8/";
	protected static String parameterFile = outputDirectory + "Parameters.txt";

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
	protected static float oxygenBulkConcentration; // = 8e1f; // [gO/um^3]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h 8.33e6
	private static float oxygenDiffusivity;// = 4e4f; // [um^2/h]

	// Toxin (T) - Toxin bulk concentration (0 initial)-V.B.
	protected static float ToxinBulkConcentration = 0.0f; // [g/um^3]
	// Toxin diffusivity-V.B. Starting with same diffusivity as oxygen
	private static float ToxinDiffusivity = 0.03f; // [um^2/h]
	private static float Toxinthreshold;

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 1.5e-7f; // [g/um^3]
	// protected static float specificMassX = 2.0e-7f; // [g/um^3]

	// Yield coefficients
	private static float Y = 0.5f; // [gX/gX]
	private static float YO = 0.5f; // [gO/gX]

	// Processes
	// Growth (biomass production)
	protected static float qsMax = 2f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]

	// Toxin production
	protected static float qT = 0.1f; // [1/h]
	protected static float f = 0.1f;
	protected static float kd1 = 0.001f;
	protected static float kd2 = 0.001f;
	protected static float Deltaf= 0.1f;
	protected static float Toxin_Investment1;
	protected static float Toxin_Investment2;
	// Toxin killing rate per Toxin unit concentration
	protected static float kT = 5;// [um^3/gT/hr]
	// Toxin loss rate per unit of susceptible species concentration
	protected static float BetaT = 0.0f; // [um^3/ gX /hr]
	protected static float alpha = 1.0f; // (gT/gX) Toxin to Biomass conversion
											// factor
	private static float KO = 3.5e-14f; // [gO/um^3]
	// Computation parameters
	protected static float systemSize = 250f; // [um]
	protected static float CompulsoryTime;
	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 1f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 25f / systemSize;

	// other model parameters
	protected static int gridSide = 65; // mult0igrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0; // NO DETACHMENT PRESENT IN THIS CASE
	protected static float h_star; //[um]
	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 130;
	protected static float InitialRatio;

	// max cell mass
	protected static float cellmassmax;

	// min cell mass
	protected static float cellmassmin;

	// biofilm carachteristic length
	protected static float Lbac;

	// delta active layer depth
	protected static double Delta;

	private static BiomassSpecies Producer1;
	private static BiomassSpecies Producer2;
	private static boolean nutrientOff; 

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
		SoluteSpecies Toxin1 = new SoluteSpecies("Toxin1", ToxinDiffusivity);
		// set up the simplest type of bulk concentration: constant
		Toxin1.setBulkConcentration(new ConstantBulkConcentration(
				ToxinBulkConcentration));
		
		// Toxin species
		SoluteSpecies Toxin2 = new SoluteSpecies("Toxin2", ToxinDiffusivity);
		// set up the simplest type of bulk concentration: constant
		Toxin2.setBulkConcentration(new ConstantBulkConcentration(
				ToxinBulkConcentration));

		// 2. Create the particulate species (soliids)
		// X active mass
		ParticulateSpecies[] activeX = new ParticulateSpecies[2];
		activeX[0] = new ParticulateSpecies("activeX",specificMassX, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		float[] fractionalVolumeCompositionH1 = new float[2];
		fractionalVolumeCompositionH1[0] = 1.0f;

		// 2.1 V.B. Create the particulate species (soliids)
		// X2 active mass
		ParticulateSpecies[] activeX2 = new ParticulateSpecies[2];
		activeX2[0] = new ParticulateSpecies("activeX2",specificMassX, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		float[] fractionalVolumeCompositionH2 = new float[2];
		fractionalVolumeCompositionH2[0] = 1.0f;

		
		ParticulateSpecies inertX = new ParticulateSpecies("inertX",specificMassX, Color.gray);
		activeX[1] = inertX;
		ParticulateSpecies inertX2 = new ParticulateSpecies("inertX2",specificMassX, Color.gray);
		activeX2[1] = inertX2;
		
		
		// 3. Create the biomass species
		// ToxinMinus = new BiomassSpecies("ToxinMinus", spX,

		// fractionalVolumeCompositionH1);

		// -----change color from Blue to Green if Toxin concentration is the
		// above threshold----

		Toxin_Investment1 = f * qsMax * Y;
		Toxin_Investment2 = (f+Deltaf) * qsMax * Y;
		ProcessFactor mO = new Saturation(oxygen, KO);

		Producer1 = new BiomassSpecies("Producer1", activeX,
				fractionalVolumeCompositionH1);
		Producer1.setActiveMass(activeX[0]);
		Producer1.setEpsMass(inertX);

		// 3.1 V.B. Create the second biomass specie
		Producer2 = new BiomassSpecies("Producer2", activeX2,
				fractionalVolumeCompositionH2);
		Producer2.setActiveMass(activeX2[0]);
		Producer2.setEpsMass(inertX2);
		
		

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growth;
		if (nutrientOff)
			growth = new Reaction("growth", activeX[0], qsMax, 0);
		else {
			growth = new Reaction("growth", activeX[0], qsMax, 1);
			growth.addFactor(mO);
		}
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		// similarly for activeX2
		Reaction growthX2;
		if (nutrientOff)
			growthX2 = new Reaction("growthX2", activeX2[0], qsMax, 0);
		else {
			growthX2 = new Reaction("growthX2", activeX2[0], qsMax, 1);
			growthX2.addFactor(mO);
		}	
		// Toxin production
		// Reaction ToxinProduction = new Reaction("ToxinProduction", activeX2,
		// qT, 1);
		Reaction ToxinProduction1;
		if (nutrientOff)
			ToxinProduction1 = new Reaction("ToxinProduction1", activeX[0],
					Toxin_Investment1, 1);
		else {
				ToxinProduction1 = new Reaction("ToxinProduction1", activeX[0],
				Toxin_Investment1, 1);
				ToxinProduction1.addFactor(mO);
		}		
		
		Reaction ToxinProduction2;
		if (nutrientOff)
			ToxinProduction2 = new Reaction("ToxinProduction2", activeX2[0],
					Toxin_Investment2, 1);
		else {
				ToxinProduction2 = new Reaction("ToxinProduction2", activeX2[0],
				Toxin_Investment2, 1);
				ToxinProduction2.addFactor(mO);
		}	
		
		// Toxin Killing effect
		Reaction ToxinKilling1 = new Reaction("ToxinKilling1", activeX[0], kT, 1);
		ProcessFactor T1 = new Linear1(Toxin2);
		ToxinKilling1.addFactor(T1);

		// Toxin Killing effect
		Reaction ToxinKilling2 = new Reaction("ToxinKilling2", activeX2[0], kT, 1);
		ProcessFactor T2 = new Linear1(Toxin1);
		ToxinKilling2.addFactor(T2);
		
		Reaction BiomassDecay1 = new Reaction ("BiomassDecay1", activeX[0], kd1, 1);
		Reaction BiomassDecay2 = new Reaction ("BiomassDecay2", activeX2[0], kd2, 1);
		
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(4);
		rsXactive.addReaction(growth, Y);
		rsXactive.addReaction(ToxinProduction1, -1);
		rsXactive.addReaction(BiomassDecay1, -1);
		rsXactive.addReaction(ToxinKilling1, -1);
		activeX[0].setProcesses(rsXactive);

		// similarly for activeX2 growth and investement in Toxin production
		NetReaction rsXactiveX2 = new NetReaction(4);
		rsXactiveX2.addReaction(growthX2, Y);
		rsXactiveX2.addReaction(BiomassDecay2, -1);
		rsXactiveX2.addReaction(ToxinProduction2, -1);
		rsXactiveX2.addReaction(ToxinKilling2, -1);
		activeX2[0].setProcesses(rsXactiveX2);
			
		// inert
		NetReaction rsXinert1 = new NetReaction(2);
		rsXinert1.addReaction(ToxinKilling1, 1);
		rsXinert1.addReaction(BiomassDecay1, 1);
		inertX.setProcesses(rsXinert1);

		// inert
		NetReaction rsXinert2 = new NetReaction(2);
		rsXinert2.addReaction(ToxinKilling2, 1);
		rsXinert2.addReaction(BiomassDecay2, 1);
		inertX2.setProcesses(rsXinert2);


		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(growth, -(YO * (1 - Y)));
		rsOxygen.addReaction(growthX2, -(YO * (1 - Y)));
		oxygen.setProcesses(rsOxygen);

		// Toxin
		NetReaction rsToxin1= new NetReaction(1);
		rsToxin1.addReaction(ToxinProduction1, alpha);
		Toxin1.setProcesses(rsToxin1);

		// Toxin
		NetReaction rsToxin2= new NetReaction(1);
		rsToxin2.addReaction(ToxinProduction2, alpha);
		Toxin2.setProcesses(rsToxin2);

		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(Producer1);
		addBiomassSpecies(Producer2);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(Toxin1);
		addSoluteSpecies(Toxin2);
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
		int[] nCells = { (int) ((1 - InitialRatio) * initialParticleNumber),
				(int) (InitialRatio * initialParticleNumber) };
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
		if (args.length < 17) {
			throw new RuntimeException(
					"input arguments missing: \n"
							+ "1: output directory (CAUTION!!! directory will be erased \n"
							+ "2: seed for random number generator \n"
							+ "3: flag for running with graphics (1 on, 0 off) \n"
							+ "4: Kt \n" 
							+ "5: Toxin Diffusivity \n"
							+ "6: 2D/3D \n"  
							+ "7: Toxin investment f (Lower) \n"
							+ "8: Oxygen Bulk Concentration \n"
							+ "9: Oxygen Diffusivity \n"
							+ "10: Initial Fraction producers \n"
							+ "11: alpha \n" 
							+ "12: Biomass density \n"
							+ "13: nutrient switch (1 on, 0 off) \n"
							+ "14: Deltaf \n" 
							+ "15: h_star \n" 
							+ "16: Kd"
							+ "17: Compulsory Time \n");
		}
		// parse inputs
		outputDirectory = args[0];
		parameterFile = outputDirectory + "Parameters.txt";
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		kT = Float.parseFloat(args[3]);
		ToxinDiffusivity = Float.parseFloat(args[4]);
		geometry = Integer.parseInt(args[5]);
		f = Float.parseFloat(args[6]);
		oxygenBulkConcentration = Float.parseFloat(args[7]);
		oxygenDiffusivity = Float.parseFloat(args[8]);
		InitialRatio = Float.parseFloat(args[9]);
		alpha = Float.parseFloat(args[10]);
		specificMassX = Float.parseFloat(args[11]);
		nutrientOff = Boolean.parseBoolean(args[12]);
		Deltaf = Float.parseFloat(args[13]);
		h_star = Float.parseFloat(args[14]);
		kd1=  Float.parseFloat(args[15]);
		CompulsoryTime=  Float.parseFloat(args[16]);
		kd2=kd1;
		rdetach = Y*qsMax*specificMassX / h_star; // 0 = NO DETACHMENT PRESENT IN THIS CASE
		
		
		cellmassmax = 3.1415f * (relativeMaximumRadius * systemSize)
				* (relativeMaximumRadius * systemSize) * systemSize / gridSide
				* specificMassX;

		cellmassmin = cellmassmax / 2.0f;

		Lbac = (alpha * kT * (cellmassmax + cellmassmin) / 2.0f)
				/ (2.0f * 3.1415f * ToxinDiffusivity);

		Toxin_Investment1 = f * qsMax * Y;
		Toxin_Investment2 = (f+Deltaf) * qsMax * Y;

		Float delta_square = new Float(oxygenBulkConcentration
				* oxygenDiffusivity * Y / (YO * (1 - Y)))
				/ ((qsMax * Y) * specificMassX
						* (relativeBoundaryLayer * systemSize) * (relativeBoundaryLayer * systemSize));
		Delta = Math.sqrt(delta_square.doubleValue());
		//
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Example9ToxinDiversity();
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
		System.out.println("Delta is:" + Delta);
		System.out.println("Toxin Investment1 is:" + Toxin_Investment1);
		System.out.println("Toxin Investment2 is:" + Toxin_Investment2);
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
			Parameters.writeToFile("Toxin investment1: " + Toxin_Investment1);
			Parameters.writeToFile("Toxin investment2: " + Toxin_Investment2);
			Parameters.writeToFile("Delta: " + Delta);
			// /
			// app.addTimedStateWriter(new PovRayWriter());
			 app.addTimedStateWriter(new SoluteConcentrationWriter());
			 app.addTimedStateWriter(new SolidsConcentrationWriter());
			 app.addTimedStateWriter(new ParticlePositionWriter());
			// app.addStateWritter(new DetachmentLevelSetWriter());
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
			spw.addSeries(new InverseOfDistancesSeries(Producer1, Producer2));
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
			//Model.model().setMaximumTimeStep(1);
			// Model.model().setFinishIterationTime(600);
			//Model.model().setMaxHeight(150);
			//Model.model().setCompulsoryTimeStep(Float.POSITIVE_INFINITY);
			Model.model().setCompulsoryTimeStep(CompulsoryTime);
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
		System.out.println("Writing the output at end of simulation");
		try {
			app.forceWriteTimedStateWriters();
		} catch (ModelException e) {
			e.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}