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
public class Example11MultiProducersToxinDiversity extends ModelHandler {
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
	protected static int numberofCellGroups = 2;

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
	/**
	 * @throws ModelException
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// Oxygen species
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		
		// Toxin species
		SoluteSpecies[] Toxin = new SoluteSpecies[numberofCellGroups];
		for (int i = 0; i < Toxin.length; i++) {
			Toxin[i] = new SoluteSpecies("Toxin" + i, ToxinDiffusivity);
			Toxin[i].setBulkConcentration(new ConstantBulkConcentration(
					ToxinBulkConcentration));
		}
   
		// 2. Create the particulate species (soliids)
		// X active mass		
		ParticulateSpecies[] activeX = new ParticulateSpecies[numberofCellGroups+1];
		for (int i = 0; i < numberofCellGroups; i++) {
			float f = ((float)(i)/(float)(numberofCellGroups));
			Color c = ColorMaps.getFullJetColor(f);
			activeX[i] =  new ParticulateSpecies("activeX"+i,
					specificMassX, c);
		}
		
		// Removed the Sensitives(only multiple producers competing)
		//activeX[numberofCellGroups-1]=new ParticulateSpecies("sensitive",
		//		specificMassX, Color.black);
		
		// inert mass 
		ParticulateSpecies inertX = new ParticulateSpecies("inertX",specificMassX, Color.gray);
		activeX[numberofCellGroups] = inertX;
		
		// 3. Create the biomass species
		ProcessFactor mO = new Saturation(oxygen, KO);
		
		BiomassSpecies[] specieArray= new BiomassSpecies[numberofCellGroups];  
		for (int i = 0; i < numberofCellGroups; i++) {
			ParticulateSpecies[] speciesArray = {activeX[i], activeX[numberofCellGroups]};
			float[] fractionArray = {1f, 0f};
			specieArray[i] = new BiomassSpecies("Species"+ i, speciesArray,
					fractionArray);
			specieArray[i].setActiveMass(activeX[i]);
			specieArray[i].setEpsMass(inertX);
		}

		// 4. Create the Reaction factors, Monod and inhibition coefficients;
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
				
		Reaction[] growth = new Reaction[numberofCellGroups];
		for (int i = 0; i < numberofCellGroups; i++) {
			if (nutrientOff) {
				growth[i] = new Reaction("growth"+i, activeX[i], qsMax, 0);
			}else{
				growth[i] = new Reaction("growth"+i, activeX[i], qsMax, 1);
				growth[i].addFactor(mO);
			}
		}

		// Toxin production
		// Reaction ToxinProduction = new Reaction("ToxinProduction", activeX2,
		// qT, 1);
		Reaction[] toxinProduction = new Reaction[numberofCellGroups];
		for (int i = 0; i < numberofCellGroups; i++) {
			if (nutrientOff) {
				toxinProduction[i] = new Reaction("ToxinProduction"+i, activeX[i],
						Toxin_Investment2, 0);
			}
			else {
				toxinProduction[i] = new Reaction("ToxinProduction"+i, activeX[i],
						Toxin_Investment2, 1);
				toxinProduction[i].addFactor(mO);
			}
		}

		// Toxin Killing effect
		Reaction[] toxinKilling = new Reaction[numberofCellGroups];
		for (int i = 0; i < numberofCellGroups; i++) {
			toxinKilling[i] = new Reaction("ToxinKillingOf"+i, activeX[i], kT, numberofCellGroups);
			for (int j = 0; j < numberofCellGroups; j++) {
				if (i!=j)
					toxinKilling[i].addFactor(new Linear1(Toxin[j]));
			}
		}
		
		// Endogenous Decay
		Reaction[] BiomassDecay = new Reaction [numberofCellGroups];
		for (int i = 0; i < numberofCellGroups; i++) {
			BiomassDecay[i] = new Reaction ("BiomassDecay"+i, activeX[i], kd1, 0);
		}

		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction[] rsXactive = new NetReaction[numberofCellGroups];
		for (int i = 0; i < numberofCellGroups; i++) {
			rsXactive[i] = new NetReaction(4);
			rsXactive[i].addReaction(growth[i], Y);
			rsXactive[i].addReaction(toxinKilling[i], -1);
			rsXactive[i].addReaction(toxinProduction[i], -1);
			rsXactive[i].addReaction(BiomassDecay[i], -1);
			activeX[i].setProcesses(rsXactive[i]);
		}
	
		// inert
		NetReaction[] rsXinert = new NetReaction[numberofCellGroups]; 
		for (int i = 0; i < numberofCellGroups; i++) {
			rsXinert[i] = new NetReaction(2);
			rsXinert[i].addReaction(toxinKilling[i], 1);
			rsXinert[i].addReaction(BiomassDecay[i], 1);
			inertX.setProcesses(rsXinert[i]);
		}

		// oxygen
		NetReaction rsOxygen = new NetReaction(numberofCellGroups);
		for (int i = 0; i < numberofCellGroups; i++) {
			rsOxygen.addReaction(growth[i], -(YO * (1 - Y)));	
		}
		oxygen.setProcesses(rsOxygen);
		
		// Toxin
		NetReaction[] rsToxin = new NetReaction[numberofCellGroups]; 
		for (int i = 0; i < numberofCellGroups; i++) {     
		rsToxin[i]=new NetReaction(1);	
		rsToxin[i].addReaction(toxinProduction[i], alpha);
		Toxin[i].setProcesses(rsToxin[i]);
		}

		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addSoluteSpecies(oxygen);
		for (int i = 0; i < numberofCellGroups; i++){
			addBiomassSpecies(specieArray[i]);
			addSoluteSpecies(Toxin[i]);
		}
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
		//int[] nCells = { (int) ((1 - InitialRatio) * initialParticleNumber),
		//		(int) (InitialRatio * initialParticleNumber) };
		//inoculateRandomly(nCells);
		int[] nCells = new int[numberofCellGroups];
		for (int i = 0; i < numberofCellGroups; i++) {
			nCells[i] = 13;
		}
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
		ApplicationComponent app = new Example11MultiProducersToxinDiversity();
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
			//spw.addSeries(new InverseOfDistancesSeries(Producer1, Producer2));
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