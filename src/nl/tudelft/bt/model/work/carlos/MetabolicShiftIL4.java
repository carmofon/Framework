
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
import nl.tudelft.bt.model.particlebased.tumor.TumorModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * @author Carlos Carmona Fontaine (carmonac@mskcc.org)
 */
public class MetabolicShiftIL4 extends TumorModelHandler {
//	public class Example1KWithTumorHandler extends TumorModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/CarmonaFontaine/Documents/Xlab/Models/delete";

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

///////////// Solute parameters	// Substrate (S) - the only solute species used here
	//protected static float substrateBulkConcentration = 1e-4f; // [g/L]
	//glucose plasma concentration = 5.5mM = 0.99g/L i.e. 1g/L
	protected static float glucoseBulkConcentration = 1f; // [g/L]
	protected static float oxygenBulkConcentration = 1f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	// Glucose diffusivity (Casciari et al 1988) 1.1e-6 cm2/s = 3.96e5um2/h mouse; 5.5e-7-2.3e-7 cm2/s = 1.98e4um2/h to 8.2e3 um2/h
	private static float glucoseDiffusivity = 8.2e4f ; // [um^2/h]
	private static float oxygenDiffusivity = 8.2e3f ; // [um^2/h]
	private static float egfDiffusivity = 4.8e4f ; // [um^2/h] (bionumbers 1.34 e-6cm2/s)
	private static float csf1Diffusivity = 4.8e4f ; // [um^2/h]
	private static float IL4Diffusivity = 1.7e4f ; // [um^2/h] ref. 4e-7cm2/s ("Diffusion of interleukin-2 from cells overlaid with cytocompatible enzyme-crosslinked gelatin hydrogels.")

	// 
///////////// Cell parameters
	

	protected static int numberofCellGroups = 2;
	//
	// Particulate species (biomass X)
	protected static float specificMassX1 = 100f; // [g/L]
	protected static float specificMassX2 = 100f; // [g/L]
	
	// Death rates (maintenance cost)
	protected static float kd1 = 0.01f;
	protected static float kd2 = 0.01f;

	// Yield coefficients ?????
	private static float oYXS = 0.045f; // [gX/gS]
	private static float Oshift = oxygenBulkConcentration/2;// K value of the switch between aerobic and anaerobic metabolisms

	// Processes
	// Growth (biomass production)
	protected static float uMax1 = 0.1f; // [1/h]
	protected static float uMax2 = 0.1f; //[1/h]

	private static float gKS = 3.5e-4f; // [g/L]
//	private static float oKS = 3.5e-4f; // [g/L]
	
	// Growth factor effect constants:
	protected static float Kegf = 0.5f;
	protected static float Kcsf1 = 0.5f;
	protected static float KIL4 = 0.001f;
	
	protected static float anaEff = 1f;// ratio between energetic rate with aerobic versus anaerobic metabolism (it changes th MuMax)

	//growth factor investment
	protected static float iegf = 0.1f;
	protected static float icsf1 = 0.1f;
	protected static float iIL4 = 0.1f;

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
	protected static float relativeBoundaryLayer = 50f/systemSize;

	// other model parameters
	protected static int gridSide = 129; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 1e-5f; // 

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 100;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// substrate
		SoluteSpecies glucose = new SoluteSpecies("glucose",
				glucoseDiffusivity);
		// set up the simplest type of bulk concentration: constant
		glucose.setBulkConcentration(new ConstantBulkConcentration(
				glucoseBulkConcentration));
		
		
		SoluteSpecies oxygen = new SoluteSpecies("oxygen",
				oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// EGF
		SoluteSpecies egf = new SoluteSpecies("egf", egfDiffusivity);
			egf.setBulkConcentration(new ConstantBulkConcentration(0));// it is only produced
		// CSF-1
		SoluteSpecies csf1 = new SoluteSpecies("csf1", csf1Diffusivity);
			csf1.setBulkConcentration(new ConstantBulkConcentration(0));
		// IL4
			SoluteSpecies IL4 = new SoluteSpecies("IL4", IL4Diffusivity);
				IL4.setBulkConcentration(new ConstantBulkConcentration(0));

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
		//speciesX2.getColorFromGrowth();
	
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		
		//	Growth factor production
		
		Reaction egfsecretion = new Reaction("egfsecretion", activeX2, iegf, 0); //species or active?
		Reaction csf1secretion = new Reaction("csf1secretion", activeX1, icsf1, 0);
		Reaction IL4secretion = new Reaction("IL4secretion", activeX1, iIL4, 0);

		//ProcessFactor mSo = new Saturation(glucose, oKS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//growth
		ProcessFactor mSg = new Saturation(glucose, gKS);
		//aerobic growth weight
		ProcessFactor aeroWeight= new Saturation(oxygen,Oshift);
		//anaerobic growth weight
		ProcessFactor anaWeight= new Inhibition(oxygen,Oshift);
		//inhibition of anaerobic growth due IL4
		ProcessFactor IL4effect= new Inhibition(IL4,KIL4);

		
		// 5. Create the reactions
		// Aerobic growth
		Reaction aeroGrowthX1 = new Reaction("aeroGrowthX1", activeX1, uMax1, 2);
		aeroGrowthX1.addFactor(mSg);
		aeroGrowthX1.addFactor(aeroWeight);

		
		Reaction aeroGrowthX2 = new Reaction("aeroGrowthX2", activeX2, uMax2, 2);
		aeroGrowthX2.addFactor(mSg);
		aeroGrowthX2.addFactor(aeroWeight);

		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		
		//Anaerobic growth
		Reaction anaGrowthX1 = new Reaction("anaGrowthX1", activeX1, anaEff * uMax1, 2);
		anaGrowthX1.addFactor(mSg);
		anaGrowthX1.addFactor(anaWeight);
		
		Reaction anaGrowthX2 = new Reaction("anaGrowthX2", activeX2, anaEff * uMax2, 3);
		anaGrowthX2.addFactor(mSg);
		anaGrowthX2.addFactor(anaWeight);
		anaGrowthX2.addFactor(IL4effect);

		
		// Endogenous Decay
		
		ProcessFactor EGFeffect= new Inhibition(egf,Kegf);
			Reaction BiomassDecayX1 = new Reaction("BiomassDecayX1", activeX1, kd1, 1);		
				BiomassDecayX1.addFactor(EGFeffect);

		ProcessFactor Csf1effect= new Inhibition(csf1,Kcsf1);
			Reaction BiomassDecayX2 = new Reaction("BiomassDecayX2", activeX2, kd2, 1);
				BiomassDecayX2.addFactor(Csf1effect);


		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactiveX1 = new NetReaction(3);
		rsXactiveX1.addReaction(aeroGrowthX1, 1);
		rsXactiveX1.addReaction(anaGrowthX1, 1);
		rsXactiveX1.addReaction(BiomassDecayX1, -1);
		activeX1.setProcesses(rsXactiveX1);
		
		NetReaction rsXactiveX2 = new NetReaction(3);
		rsXactiveX2.addReaction(aeroGrowthX2, 1);
		rsXactiveX2.addReaction(anaGrowthX2, 1);
		rsXactiveX2.addReaction(BiomassDecayX2, -1);
		activeX2.setProcesses(rsXactiveX2);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsGlucose = new NetReaction(4);
		rsGlucose.addReaction(aeroGrowthX1, -(1 / oYXS));
		rsGlucose.addReaction(aeroGrowthX2, -(1 / oYXS));
		rsGlucose.addReaction(anaGrowthX1, -(1 / oYXS));
		rsGlucose.addReaction(anaGrowthX2, -(1 / oYXS));
		glucose.setProcesses(rsGlucose);
		
		NetReaction rsOxygen = new NetReaction(4);
		rsOxygen.addReaction(aeroGrowthX1, -(1 / oYXS));
		rsOxygen.addReaction(aeroGrowthX2, -(1 / oYXS));
		rsOxygen.addReaction(anaGrowthX1, -(1 / oYXS));
		rsOxygen.addReaction(anaGrowthX2, -(1 / oYXS));
		oxygen.setProcesses(rsOxygen);
		
		NetReaction rsEgf = new NetReaction(3);
		rsEgf.addReaction(egfsecretion, 1);		
		rsEgf.addReaction(aeroGrowthX1, -(1 / oYXS));
		rsEgf.addReaction(anaGrowthX1, -(1 / oYXS));
		egf.setProcesses(rsEgf);
		
		NetReaction rsCsf1 = new NetReaction(3);
		rsCsf1.addReaction(csf1secretion, 1);		
		rsCsf1.addReaction(aeroGrowthX2, -(1 / oYXS));
		rsCsf1.addReaction(anaGrowthX2, -(1 / oYXS));
		csf1.setProcesses(rsCsf1);
		
		NetReaction rsIL4 = new NetReaction(3);
		rsIL4.addReaction(IL4secretion, 1);		
		rsIL4.addReaction(aeroGrowthX2, -(1 / oYXS));
		rsIL4.addReaction(anaGrowthX2, -(1 / oYXS));
		IL4.setProcesses(rsIL4);
		
		
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesX1);
		addBiomassSpecies(speciesX2);
		addSoluteSpecies(glucose);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(egf);
		addSoluteSpecies(csf1);
		addSoluteSpecies(IL4);

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
		ApplicationComponent app = new MetabolicShiftIL4();
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
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter());
			//	app.addTimedStateWriter(new ImageWriter(glucose));
			//	app.addTimedStateWriter(new ImageWriter(oxygen));
			//	app.addTimedStateWriter(new ImageWriter(egf));
			//	app.addTimedStateWriter(new ImageWriter(csf1));
			//	app.addTimedStateWriter(new ImageWriter(IL4));
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