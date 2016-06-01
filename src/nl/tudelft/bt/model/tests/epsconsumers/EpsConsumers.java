package nl.tudelft.bt.model.tests.epsconsumers;

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
 * Test model for eps producers and the competition between bacteria that
 * produce eps hydrolysis enzymes (altruists) and egoist enzymes that produce
 * and consume EPS but do not produce it
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class EpsConsumers extends ModelHandler {
	//output directory name
	private static String outputDirectory = "C:/joao/results/lixo";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "E:/results/epsConsumers/" directory
	// "E:/results/" must exist in your computer.
	// The directory "E:/results/epsConsumers/" will be created
	// if it is non-existant. If it exists, its contents will be deleted
	// during the program initialization

	// geometry of simulation (2D or 3D)
	private static int geometry = 2;

	//Parameters for solute species
	//Substrate
	private static float substrateBulkConcentration = 0.1f; //[gCOD-S/l]

	private static float substrateSourceDiffusivity = 4e6f; //[um2/h]

	private static float KS = 4e-3f; //[gO/l]

	// duration of feast and famine cycles
	private static float feastTime = 4f; //[h]

	private static float famineTime = 24f - feastTime; //[h]

	private static float initialFeastTime = feastTime; //[h]

	//Oxygen
	private static float oxygenBulkConcentration = 4e-3f; //[gO/l]

	private static float oxygenDiffusivity = 8e6f; //[um2/h]

	private static float KO = 3.5e-4f; //[gO/l]

	private static float KO_star = 7E-04f; //[gO/l]

	//EPS hydrolysis enzyme (Enz)
	private static float enzBulkConcentration = 0f; //[gCOD-Enz/L]

	private static float enzDiffusivity = 8e4f; //[um2/h]

	//Parameters for particulate species
	//density of active biomass (either H-PHB or H-EPS)
	private static float densityH = 200.0f; //[gCOD-X/l]

	//density of EPS
	private static float densityEPS = densityH / 6; //[gCOD-EPS/l]

	//Yield coefficients
	//biomass on substrate
	private static float YXS = 0.495f; //[gCOD-X/gCOS-S]

	//polymer on substrate consumed
	private static float YPS = 0.667f; //[gCOD-EPS/gCOS-S]

	//enzyme on active mass consumed
	private static float YEnzX = 5f; //[gCOD-X/gCOS-PHB]

	// Processes
	private static float qMax = 0.952f; //[gCOD-S/gCOD-X/h]

	private static float qEnz = 1f; //[gCOD-Enz/gCOD-X/h]

	private static float uMax = YXS * qMax; //[gCOD-X/gCOD-X/h]

	// Polymer inhibition
	private static float KEPS = 1e-4f; //[fP]

	// Computation parameters
	private static float systemSize = 1000; // [micron]

	// maximum radius of biomass particles, relative to the system size
	private static float relativeMaximumRadius = 0.007f;

	// minimum radius of biomass particles, relative to the system size
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// boundary layer thickness, relative to the system size
	private static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	private static int gridSide = 33; // multigrid grid side

	private static float kShov = 1.0f; // shoving parameter

	// detachment rate coefficient
	private static float rdetach = 3e-3f;// detachment constant

	//inoculation
	private static int initialCellNumber = 31;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		//           ---- Define the solutes ----
		//---substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateSourceDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		//The following code is used for feast/famine cycles to
		//keep acetate concentration average in time. To use, comment the
		//lines above and uncomment the lines below.
		//BEGIN
		//float ac = substrateBulkConcentration
		//				/ (feastTime / (feastTime + famineTime));
		//substrate.setBulkConcentration(new ItermittentBulkConcentration(ac,
		//				initialFeastTime, feastTime, famineTime));
		//END
		//---oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		//---EPS hydrolysis enzyme
		SoluteSpecies enz = new SoluteSpecies("Enzyme", enzDiffusivity);
		enz.setBulkConcentration(new ConstantBulkConcentration(
				enzBulkConcentration));
		//           ---- Create the particulates ----
		//---X1 active mass
		ParticulateSpecies activeX1 = new ParticulateSpecies("X1", densityH,
				Color.blue);
		//---X2 active mass
		ParticulateSpecies activeX2 = new ParticulateSpecies("X2", densityH,
				Color.red);
		//---EPS-X1
		ParticulateSpecies epsX1 = new ParticulateSpecies("epsX1", densityEPS,
				Color.yellow);
		//---EPS-X2
		ParticulateSpecies epsX2 = new ParticulateSpecies("epsX2", densityEPS,
				Color.yellow);
		//           ---- Create the biomass species ----
		//----PHB producer
		//array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spX1 = { activeX1, epsX1 };
		float[] fractionalVolumeCompositionX1 = { 1.0f, 0 };
		BiomassSpecies speciesX1 = new BiomassSpecies("X1", spX1,
				fractionalVolumeCompositionX1);
		speciesX1.setActiveMass(activeX1);
		speciesX1.setEpsMass(epsX1);
		//----EPS producer
		//array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spX2 = { activeX2, epsX2 };
		float[] fractionalVolumeCompositionX2 = { 1.0f, 0 };
		BiomassSpecies speciesX2 = new BiomassSpecies("X2", spX2,
				fractionalVolumeCompositionX2);
		speciesX2.setActiveMass(activeX2);
		speciesX2.setEpsMass(epsX2); //defines the species that is modeled as
		// capsule
		//           ---- Create the process terms ----
		//Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor iS = new Inhibition(substrate, KS);
		ProcessFactor mO = new Saturation(oxygen, KO);
		ProcessFactor mO_star = new Saturation(oxygen, KO_star);
		ProcessFactor mEPSX1 = new Saturation(epsX1, KEPS);
		ProcessFactor mEPSX2 = new Saturation(epsX2, KEPS);
		//           ---- Create the reactions ----
		//---carbon source uptake by X1
		Reaction sUptakeX1 = new Reaction("substrate uptake by X1", activeX1,
				qMax, 2);
		sUptakeX1.addFactor(mS);
		sUptakeX1.addFactor(mO);
		//---acetate uptake by X2
		Reaction sUptakeX2 = new Reaction("substrate uptake by X2", activeX2,
				qMax, 2);
		sUptakeX2.addFactor(mS);
		sUptakeX2.addFactor(mO);
		//---Growth of X1
		Reaction growthX1 = new Reaction("growth X1", activeX1, uMax, 2);
		growthX1.addFactor(mS);
		growthX1.addFactor(mO_star);
		//---Growth of X2
		Reaction growthX2 = new Reaction("growth X2", activeX2, uMax, 2);
		growthX2.addFactor(mS);
		growthX2.addFactor(mO_star);
		//---Enzyme production by X1
		Reaction enzymeProduction = new Reaction("enzymeProduction", activeX1,
				qEnz, 2);
		enzymeProduction.addFactor(mS);
		enzymeProduction.addFactor(mO_star);
		// 		---- Assign reaction to each species through NetRaction ----
		//---X1
		NetReaction rsX1 = new NetReaction(2); 
		rsX1.addReaction(growthX1, 1);
		rsX1.addReaction(enzymeProduction, -YEnzX);
		activeX1.setProcesses(rsX1);
		//---X2 active mass
		NetReaction rsX2 = new NetReaction(1);
		rsX2.addReaction(growthX2, 1);
		activeX2.setProcesses(rsX2);
		//---EPS X1
		NetReaction rsEpsX1 = new NetReaction(2);
		rsEpsX1.addReaction(sUptakeX1, YPS);
		rsEpsX1.addReaction(growthX1, -YPS / YXS);
		epsX1.setProcesses(rsEpsX1);
		//---EPS X2
		NetReaction rsEpsX2 = new NetReaction(2);
		rsEpsX2.addReaction(sUptakeX2, YPS);
		rsEpsX2.addReaction(growthX2, -YPS / YXS);
		epsX2.setProcesses(rsEpsX2);
		//---Substrate
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(sUptakeX1, -1);
		rsSubstrate.addReaction(sUptakeX2, -1);
		substrate.setProcesses(rsSubstrate);
		//---oxygen
		NetReaction rsOxygen = new NetReaction(5);
		rsOxygen.addReaction(sUptakeX1, -(1 - YPS));
		rsOxygen.addReaction(growthX1, -(YPS / YXS - 1));
		rsOxygen.addReaction(sUptakeX2, -(1 - YPS));
		rsOxygen.addReaction(growthX2, -(YPS / YXS - 1));
		rsOxygen.addReaction(enzymeProduction, -(YEnzX - 1));
		oxygen.setProcesses(rsOxygen);
		//---Enzyme
		NetReaction rsEnzyme = new NetReaction(1);
		rsEnzyme.addReaction(enzymeProduction, 1);
		enz.setProcesses(rsEnzyme);
		//           ---- Add the species to system ----
		addBiomassSpecies(speciesX1);
		addBiomassSpecies(speciesX2);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(substrate);
		addSoluteSpecies(enz);
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
		try {
			int[] nCells = { initialCellNumber, initialCellNumber };
			inoculateRandomly(nCells);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
			return;
		}
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// set preprocessing and postprocessing steps for multigrid algorithm
		MultigridVariable.setSteps(10, 100);
		//
		ApplicationComponent app = new EpsConsumers();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		//uncomment the following line for plot of produced biomass
		//app = new SeriesVizualizer(app, prod);
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		//uncomment the following line for plot of total biomass in biofilm
		//app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		app = new DetachedBiomassVizualizer(app);
		// bulk concentrations
		app = new BulkConcentrationVizualizer(app);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			// set the global mass balance parameters
			// --- nothing to set in this case: constant bulk concentration
			//initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
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
			// write the header in parameters.txt file
			spw.initializeParametersWriting();
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			//wait for user to press start iteration
			//app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}