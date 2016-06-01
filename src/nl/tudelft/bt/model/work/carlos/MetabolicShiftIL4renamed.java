
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
public class MetabolicShiftIL4renamed extends TumorModelHandler {
//	public class Example1KWithTumorHandler extends TumorModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/Karlo11/Documents/Xlab/Model/MetabolicShift/";

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
	private static float glucoseDiffusivity = 8.2e6f ; // [um^2/h]
	private static float oxygenDiffusivity = 8.2e5f ; // [um^2/h]
	private static float egfDiffusivity = 4.8e1f ; // [um^2/h] (bionumbers 1.34 e-6cm2/s)
	private static float csf1Diffusivity = 4.8e1f ; // [um^2/h]
	private static float IL4Diffusivity = 1.7e4f ; // [um^2/h] ref. 4e-7cm2/s ("Diffusion of interleukin-2 from cells overlaid with cytocompatible enzyme-crosslinked gelatin hydrogels.")

	// 
///////////// Cell parameters
	

	protected static int numberofCellGroups = 2;
	//
	// Particulate species (biomass X)
	protected static float specificMassCC = 100f; // [g/L]
	protected static float specificMassTAM = 100f; // [g/L]
	
	// Death rates (maintenance cost)
	protected static float kdCC = 0.01f;
	protected static float kdTAM = 0.01f;

	// Yield coefficients ?????
	private static float oYXS = 0.045f; // [gX/gS]
	private static float Oshift = oxygenBulkConcentration/10;// K value of the switch between aerobic and anaerobic metabolisms

	// Processes
	// Growth (biomass production)
	protected static float uMaxCC = 0.1f; // [1/h]
	protected static float uMaxTAM= 0.1f; //[1/h]

	private static float gKS = 3.5e-4f; // [g/L]
	
//	private static float oKS = 3.5e-4f; // [g/L]
	
	// Growth factor effect constants:
	protected static float Kegf = 0.5f;
	protected static float Kcsf1 = 0.5f;
	protected static float KIL4 = 0.001f;
	
	protected static float anaEff = 1f;// ratio between energetic rate with aerobic versus anaerobic metabolism (it changes th MuMax)

	//growth factor investment
	protected static float iegf = 0.01f;
	protected static float icsf1 = 0.01f;
	protected static float iIL4 = 0.01f;

	// Computation parameters 
	protected static float systemSize = 2000; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 5f/systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.6f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 100f/systemSize;

	// other model parameters
	protected static int gridSide = 129; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 1e-5f; // 

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 250;
	
	protected static SoluteSpecies glucose;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// substrate
		glucose = new SoluteSpecies("glucose",
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
		ParticulateSpecies activeCC = new ParticulateSpecies("activeCC",
				specificMassCC, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spCC = { activeCC };
		float[] fractionalVolumeCompositionCC = { 1.0f };
		
		ParticulateSpecies activeTAM = new ParticulateSpecies("activeTAM",
				specificMassTAM, Color.green);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spTAM = { activeTAM };
		float[] fractionalVolumeCompositionTAM = { 1.0f };

		// 3. Create the biomass species
		BiomassSpecies speciesCC = new BiomassSpecies("speciesCC", spCC,
				fractionalVolumeCompositionCC);
		speciesCC.setActiveMass(activeCC);
		//speciesCC.getColorFromGrowth();

		BiomassSpecies speciesTAM = new BiomassSpecies("speciesTAM", spTAM,
				fractionalVolumeCompositionTAM);
		speciesTAM.setActiveMass(activeTAM);
		//speciesTAM.getColorFromGrowth();
	
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		
		//	Growth factor production
		
		Reaction egfsecretion = new Reaction("egfsecretion", activeTAM, iegf, 0); //species or active?
		Reaction csf1secretion = new Reaction("csf1secretion", activeCC, icsf1, 0);
		Reaction IL4secretion = new Reaction("IL4secretion", activeCC, iIL4, 0);

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
		Reaction aeroGrowthCC = new Reaction("aeroGrowthCC", activeCC, uMaxCC, 2);
		aeroGrowthCC.addFactor(mSg);
		aeroGrowthCC.addFactor(aeroWeight);

		
		Reaction aeroGrowthTAM = new Reaction("aeroGrowthTAM", activeTAM, uMaxTAM, 2);
		aeroGrowthTAM.addFactor(mSg);
		aeroGrowthTAM.addFactor(aeroWeight);

		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		
		//Anaerobic growth
		Reaction anaGrowthCC = new Reaction("anaGrowthCC", activeCC, anaEff * uMaxCC, 2);
		anaGrowthCC.addFactor(mSg);
		anaGrowthCC.addFactor(anaWeight);
		
		Reaction anaGrowthTAM = new Reaction("anaGrowthTAM", activeTAM, anaEff * uMaxTAM, 3);
		anaGrowthTAM.addFactor(mSg);
		anaGrowthTAM.addFactor(anaWeight);
		anaGrowthTAM.addFactor(IL4effect);

		
		// Endogenous Decay
		
		ProcessFactor EGFeffect= new Inhibition(egf,Kegf);
			Reaction BiomassDecayCC = new Reaction("BiomassDecayCC", activeCC, kdCC, 1);		
				BiomassDecayCC.addFactor(EGFeffect);

		ProcessFactor Csf1effect= new Inhibition(csf1,Kcsf1);
			Reaction BiomassDecayTAM = new Reaction("BiomassDecayTAM", activeTAM, kdTAM, 1);
				BiomassDecayTAM.addFactor(Csf1effect);


		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactiveCC = new NetReaction(3);
		rsXactiveCC.addReaction(aeroGrowthCC, 1);
		rsXactiveCC.addReaction(anaGrowthCC, 1);
		rsXactiveCC.addReaction(BiomassDecayCC, -1);
		activeCC.setProcesses(rsXactiveCC);
		
		NetReaction rsXactiveTAM = new NetReaction(3);
		rsXactiveTAM.addReaction(aeroGrowthTAM, 1);
		rsXactiveTAM.addReaction(anaGrowthTAM, 1);
		rsXactiveTAM.addReaction(BiomassDecayTAM, -1);
		activeTAM.setProcesses(rsXactiveTAM);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsGlucose = new NetReaction(4);
		rsGlucose.addReaction(aeroGrowthCC, -(1 / oYXS));
		rsGlucose.addReaction(aeroGrowthTAM, -(1 / oYXS));
		rsGlucose.addReaction(anaGrowthCC, -(1 / oYXS));
		rsGlucose.addReaction(anaGrowthTAM, -(1 / oYXS));
		glucose.setProcesses(rsGlucose);
		
		NetReaction rsOxygen = new NetReaction(2);
		//NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(aeroGrowthCC, -(1 / oYXS));
		rsOxygen.addReaction(aeroGrowthTAM, -(1 / oYXS));

		oxygen.setProcesses(rsOxygen);
		
		NetReaction rsEgf = new NetReaction(3);
		rsEgf.addReaction(egfsecretion, 1);		
		rsEgf.addReaction(aeroGrowthCC, -(1 / oYXS));
		rsEgf.addReaction(anaGrowthCC, -(1 / oYXS));
		egf.setProcesses(rsEgf);
		
		NetReaction rsCsf1 = new NetReaction(3);
		rsCsf1.addReaction(csf1secretion, 1);		
		rsCsf1.addReaction(aeroGrowthTAM, -(1 / oYXS));
		rsCsf1.addReaction(anaGrowthTAM, -(1 / oYXS));
		csf1.setProcesses(rsCsf1);
		
		NetReaction rsIL4 = new NetReaction(3);
		rsIL4.addReaction(IL4secretion, 1);		
		rsIL4.addReaction(aeroGrowthTAM, -(1 / oYXS));
		rsIL4.addReaction(anaGrowthTAM, -(1 / oYXS));
		IL4.setProcesses(rsIL4);
		
		
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesCC);
		addBiomassSpecies(speciesTAM);
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
		ApplicationComponent app = new MetabolicShiftIL4renamed();
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
			app.initializeDiffusionReactionSystem(); // also innoculates
			
			//grafics
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter());
			//app.addTimedStateWriter(new ImageWriter(glucose));
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