package nl.tudelft.bt.model.examples;

import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.*;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;

/**
 * The model used for the case study described in "A Framework for
 * Multidimensional Modelling of Activity and Structure of Multispecies
 * Biofilms" by Joao B. Xavier, Cristian Picioreanu and Mark C.M. van
 * Loosdrecht.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Panos extends ModelHandler {
	//output directory name
	private static String outputDirectory = "Z:/joao/scrapbook/chrysi/results";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "c:\results\phb\" directory "c:\results" must
	// exist in your computer. The directory "c:\results\phb\" will be created
	// if it is non-existant. If it exists, its contents will be deleted
	// during the program initialization

	// geometry of simulation (2D or 3D)
	private static int geometry = 2;

	//Parameters for solute species
	//Substrate
	private static float substrateBulkConcentration = 0.1f; //[gCOD_S/l]

	private static float substrateSourceDiffusivity = 4e6f; //[um2/h]

	private static float KS = 1e-3f; //[gO/l]
	private static float KU = 1e-1f; //[gO/l]
	private static float KB = 0.085f; //[gO/l]
	private static float kUap = 0.05f;
	private static float kEps = 0.18f;

	//Oxygen
	private static float oxygenBulkConcentration = 4e-3f; //[gO/l]

	private static float oxygenDiffusivity = 8e6f; //[um2/h]

	private static float KO = 3.5e-4f; //[gO/l]

	//Parameters for particulate species
	//density of active biomass (either H-PHB or H-EPS)
	private static float densityX = 200.0f; //[gCOD_X/l]

	//density of EPS
	private static float densityEPS = densityX / 2; //[gCOD_EPS/l]

	//density of inert biomass
	private static float densityI = densityX + 20; //[gCOD_I/l]

	//polymer (either PHB or EPS) on substrate
	private static float YS = 0.35f; //[gCOD_PHB/gCOS_S]
	private static float YP = 0.45f; //[gCOD_PHB/gCOS_S]

	// Processes
	//Substrate uptake
	private static float qsMax = 0.952f; //[gCOD_S/gCOD_X/h]
	private static float quMax = 0.053f; //[gCOD_S/gCOD_X/h]
	private static float qbMax = 0.003f; //[gCOD_S/gCOD_X/h]
	private static float rHydr = 0.17f / 24.0f;
	private static float fd = 0.8f;

	//Biomass decay rate constant
	private static float kDecay = 0.033f; //[gCOD_X/gCOD_X/h]

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
	private static float bDetach = 0.1f / 24.0f; // [hour^-1]

	//inoculation
	private static int initialCellNumber = 100;

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
		
		SoluteSpecies uap = new SoluteSpecies("UAP", substrateSourceDiffusivity);
		uap.setBulkConcentration(new ConstantBulkConcentration(0.0f));
		SoluteSpecies bap = new SoluteSpecies("BAP", substrateSourceDiffusivity);
		bap.setBulkConcentration(new ConstantBulkConcentration(0.0f));
		
		//           ---- Create the particulates ----
		//--- active biomass
		ParticulateSpecies active = new ParticulateSpecies(
				"active", densityX, Color.blue);
		//---EPS
		ParticulateSpecies eps = new ParticulateSpecies("eps", densityEPS,
				Color.yellow);
		//---Inert
		//Inert is descriminated by species to assess mass balances
		ParticulateSpecies inert = new ParticulateSpecies("inerts",
				densityI, Color.red);
		//           ---- Create the biomass species ----
		//----PHB producer
		//array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] heterotrophs = { active, inert, eps };
		float[] fractionalVolumeComposition = { 1.0f, 0.0f, 0.0f };
		BiomassSpecies species = new BiomassSpecies("EPS producer", heterotrophs,
				fractionalVolumeComposition);
		species.setActiveMass(active);
		species.setInertMass(inert);
		species.setEpsMass(eps);
		
		//           ---- Create the process terms ----
		//Monod coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor mO = new Saturation(oxygen, KO);
		
		ProcessFactor mU = new Saturation(uap, KU);
		ProcessFactor mB = new Saturation(bap, KB);
		//           ---- Create the reactions ----
		//---carbon source uptake by H-PHB
		
		
		Reaction sUtilization = new Reaction("substrate utilization",
				active, qsMax, 2);
		sUtilization.addFactor(mS);
		sUtilization.addFactor(mO);
		//---acetate uptake by H-EPS
		Reaction uapUtilization = new Reaction("UAP utilization",
				active, quMax, 2);
		uapUtilization.addFactor(mU);
		uapUtilization.addFactor(mO);
		
		Reaction bapUtilization = new Reaction("BAP utilization",
				active, qbMax, 2);
		bapUtilization.addFactor(mB);
		bapUtilization.addFactor(mO);
		
		Reaction epsHydrolysis = new Reaction("EPS hydrolysis",
				eps, rHydr, 0);
		
		Reaction biomassDecay = new Reaction("Biomass decay",
				active, kDecay, 0);
		
		Reaction biomassDetachment = new Reaction("Biomass detachment",
				active, bDetach, 0);
		
		// 		---- Assign reaction to each species through NetRaction ----
		//---H-PHB active mass
		
		NetReaction rs = new NetReaction(1);
		rs.addReaction(sUtilization, -1);
		substrate.setProcesses(rs);
		
		NetReaction rsOxygen = new NetReaction(3);
		rsOxygen.addReaction(sUtilization, -(1 - YS));
		rsOxygen.addReaction(uapUtilization, -(1 - YP));
		rsOxygen.addReaction(bapUtilization, -(1 - YP));
		oxygen.setProcesses(rsOxygen);
		
		NetReaction rsUAP = new NetReaction(2);
		rsUAP.addReaction(sUtilization, kUap);
		rsUAP.addReaction(uapUtilization, -1);
		uap.setProcesses(rsUAP);
		
		NetReaction rsBAP = new NetReaction(2);
		rsBAP.addReaction(bapUtilization, -1);
		rsBAP.addReaction(epsHydrolysis, 1);
		bap.setProcesses(rsBAP);
		
		NetReaction rsEPS = new NetReaction(2);
		rsEPS.addReaction(sUtilization, kEps);
		rsEPS.addReaction(epsHydrolysis, -1);
		eps.setProcesses(rsEPS);
		
		NetReaction rsActive = new NetReaction(5);
		rsActive.addReaction(sUtilization, YS * (1 - kUap - kEps));
		rsActive.addReaction(uapUtilization, YP);
		rsActive.addReaction(bapUtilization, YP);
		rsActive.addReaction(biomassDecay, -1);
		rsActive.addReaction(biomassDetachment, -bDetach);
		active.setProcesses(rsActive);
		
		NetReaction rsInert = new NetReaction(1);
		rsInert.addReaction(biomassDecay, (1 - fd));
		inert.setProcesses(rsInert);
		
		//           ---- Add the species to system ----
		addBiomassSpecies(species);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(substrate);
		addSoluteSpecies(uap);
		addSoluteSpecies(bap);
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
			int[] nCells = { initialCellNumber };
			inoculateRandomly(nCells);
		} catch (Exception e) {
			System.out.println(e);
			System.exit(-1);
			return;
		}
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new MassDetachment(0);
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
		ApplicationComponent app = new Panos();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		//uncomment the following line for plot of produced biomass
//		app = new SeriesVizualizer(app, prod);
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
			spw.addSeries(Model.model().detachedBiomassContainer()
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