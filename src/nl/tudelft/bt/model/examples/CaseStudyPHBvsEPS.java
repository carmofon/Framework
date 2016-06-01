package nl.tudelft.bt.model.examples;

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
 * The model used for the case study described in "A Framework for
 * Multidimensional Modelling of Activity and Structure of Multispecies
 * Biofilms" by Joao B. Xavier, Cristian Picioreanu and Mark C.M. van
 * Loosdrecht.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class CaseStudyPHBvsEPS extends ModelHandler {
	// output directory name
	private static String outputDirectory = "/Users/Carlos/Documents/Results/Model/Case/";

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

	// Parameters for solute species
	// Substrate
	private static float substrateBulkConcentration = 0.1f; // [gCOD_S/l]

	private static float substrateSourceDiffusivity = 4e6f; // [um2/h]

	private static float KS = 4e-3f; // [gO/l]

	// duration of feast and famine cycles
	private static float feastTime = 4f; // [h]

	private static float famineTime = 24f - feastTime; // [h]

	private static float initialFeastTime = feastTime; // [h]

	// Oxygen
	private static float oxygenBulkConcentration = 4e-3f; // [gO/l]

	private static float oxygenDiffusivity = 8e6f; // [um2/h]

	private static float KO = 3.5e-4f; // [gO/l]

	private static float KO_star = 7E-04f; // [gO/l]

	// Parameters for particulate species
	// density of active biomass (either H-PHB or H-EPS)
	private static float densityH = 200.0f; // [gCOD_X/l]

	// density of EPS
	private static float densityEPS = densityH / 6; // [gCOD_EPS/l]

	// density of PHB
	private static float densityPHB = densityH * 5; // [gCOD_PHB/l]

	// density of inert biomass
	private static float densityI = densityH; // [gCOD_I/l]

	// Yield coefficients
	// biomass on substrate
	private static float YSX = 0.495f; // [gCOD_X/gCOS_S]

	// polymer (either PHB or EPS) on substrate
	private static float YSP = 0.75f; // [gCOD_PHB/gCOS_S]

	// biomass on PHB
	private static float YPX = 0.624f; // [gCOD_X/gCOS_PHB]

	// Processes
	// Substrate uptake
	private static float qMax = 0.952f; // [gCOD_S/gCOD_X/h]

	// PHB consumption
	private static float kPhb = 0.15f; // [gCOD_PHB/gCOD_X/h]

	// Polymer inhibition
	private static float KP = 1f; // [fP]

	// Biomass decay rate constant
	private static float kDecay = 0.0033f; // [gCOD_X/gCOD_X/h]

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

	// inoculation
	private static int initialCellNumber = 31;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// ---- Define the solutes ----
		// ---substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateSourceDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		// The following code is used for feast/famine cycles to
		// keep acetate concentration average in time. To use, comment the
		// lines above and uncomment the lines below.
		// BEGIN
		// float ac = substrateBulkConcentration
		// / (feastTime / (feastTime + famineTime));
		// substrate.setBulkConcentration(new ItermittentBulkConcentration(ac,
		// initialFeastTime, feastTime, famineTime));
		// END
		// ---oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// ---- Create the particulates ----
		// ---H-PHB active mass
		ParticulateSpecies activeHPHB = new ParticulateSpecies(
				"H-PHB-active-mass", densityH, Color.blue);
		// ---H-EPS active mass
		ParticulateSpecies activeHEPS = new ParticulateSpecies(
				"H-EPS-active-mass", densityH, Color.red);
		// ---PHB
		ParticulateSpecies phb = new ParticulateSpecies("phb", densityPHB,
				Color.yellow);
		// ---EPS
		ParticulateSpecies eps = new ParticulateSpecies("eps", densityEPS,
				Color.yellow);
		// ---Inert
		// Inert is descriminated by species to assess mass balances
		ParticulateSpecies inertHPHB = new ParticulateSpecies("inert-H-PHB",
				densityI, Color.gray);
		ParticulateSpecies inertHEPS = new ParticulateSpecies("inert-H-EPS",
				densityI, Color.gray);
		// ---- Create the biomass species ----
		// ----PHB producer
		// array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spHPHB = { activeHPHB, phb, inertHPHB };
		float[] fractionalVolumeCompositionHPHB = { 1.0f, 0, 0 };
		BiomassSpecies speciesHPHB = new BiomassSpecies("PHB producer", spHPHB,
				fractionalVolumeCompositionHPHB);
		speciesHPHB.setActiveMass(activeHPHB);
		speciesHPHB.setInertMass(inertHPHB);
		// ----EPS producer
		// array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spHEPS = { activeHEPS, eps, inertHEPS };
		float[] fractionalVolumeCompositionHEPS = { 1.0f, 0, 0 };
		BiomassSpecies speciesHEPS = new BiomassSpecies("EPS producer", spHEPS,
				fractionalVolumeCompositionHEPS);
		speciesHEPS.setActiveMass(activeHEPS);
		speciesHEPS.setInertMass(inertHEPS);
		speciesHEPS.setEpsMass(eps); // defines the species that is modeled
										// as
		// capsule
		// ---- Create the process terms ----
		// Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor iS = new Inhibition(substrate, KS);
		ProcessFactor mO = new Saturation(oxygen, KO);
		ProcessFactor mO_star = new Saturation(oxygen, KO_star);
		ProcessFactor iPHB = new InhibitionFromFraction(phb, activeHPHB, KP);
		ProcessFactor iEPS = new InhibitionFromFraction(eps, activeHEPS, KP);
		// ---- Create the reactions ----
		// ---carbon source uptake by H-PHB
		Reaction sUptakeHPHB = new Reaction("substrate uptake by H-PHB",
				activeHPHB, qMax, 3);
		sUptakeHPHB.addFactor(mS);
		sUptakeHPHB.addFactor(mO);
		sUptakeHPHB.addFactor(iPHB);
		// ---acetate uptake by H-EPS
		Reaction sUptakeH2 = new Reaction("substrate uptake by H-EPS",
				activeHEPS, qMax, 3);
		sUptakeH2.addFactor(mS);
		sUptakeH2.addFactor(mO);
		sUptakeH2.addFactor(iEPS);
		// ---Growth of H-PHB
		Reaction growthHPHB = new Reaction("growth H-PHB", activeHPHB, qMax
				* YSX, 2);
		growthHPHB.addFactor(mS);
		growthHPHB.addFactor(mO_star);
		// ---Growth of H-EPS
		Reaction growthHEPS = new Reaction("growth H-EPS", activeHEPS, qMax
				* YSX, 2);
		growthHEPS.addFactor(mS);
		growthHEPS.addFactor(mO_star);
		// ---decay of H-PHB
		Reaction decayHPHB = new Reaction("decay H-PHB", activeHPHB, kDecay, 0);
		// ---decay H-EPS
		Reaction decayHEPS = new Reaction("decayH H-EPS", activeHEPS, kDecay, 0);
		// ---PHB consumption
		Reaction phbConsumption = new Reaction("PHB consumption", phb, kPhb, 2);
		phbConsumption.addFactor(iS);
		phbConsumption.addFactor(mO);
		// ---- Assign reaction to each species through NetRaction ----
		// ---H-PHB active mass
		NetReaction rsHPHBactive = new NetReaction(3); // 3 is number of
		// reactions
		rsHPHBactive.addReaction(growthHPHB, 1);
		rsHPHBactive.addReaction(phbConsumption, YPX);
		rsHPHBactive.addReaction(decayHPHB, -1);
		activeHPHB.setProcesses(rsHPHBactive);
		// ---H-EPS active mass
		NetReaction rsHEPSactive = new NetReaction(2);
		rsHEPSactive.addReaction(growthHEPS, 1);
		rsHEPSactive.addReaction(decayHEPS, -1);
		activeHEPS.setProcesses(rsHEPSactive);
		// ---PHB
		NetReaction rsPhb = new NetReaction(3);
		rsPhb.addReaction(sUptakeHPHB, YSP);
		rsPhb.addReaction(growthHPHB, -YSP / YPX);
		rsPhb.addReaction(phbConsumption, -1);
		phb.setProcesses(rsPhb);
		// ---EPS
		NetReaction rsEps = new NetReaction(2);
		rsEps.addReaction(sUptakeH2, YSP);
		rsEps.addReaction(growthHEPS, -YSP / YPX);
		eps.setProcesses(rsEps);
		// ---InertH1
		NetReaction rsIHPHB = new NetReaction(1);
		rsIHPHB.addReaction(decayHPHB, 0.4f);
		inertHPHB.setProcesses(rsIHPHB);
		// ---InertH2
		NetReaction rsIEPS = new NetReaction(1);
		rsIEPS.addReaction(decayHEPS, 0.4f);
		inertHEPS.setProcesses(rsIEPS);
		// ---carbon source
		NetReaction rsAcetate = new NetReaction(4);
		rsAcetate.addReaction(sUptakeHPHB, -1);
		rsAcetate.addReaction(sUptakeH2, -1);
		rsAcetate.addReaction(decayHPHB, +0.6f);
		rsAcetate.addReaction(decayHEPS, +0.6f);
		substrate.setProcesses(rsAcetate);
		// ---oxygen
		NetReaction rsOxygen = new NetReaction(5);
		rsOxygen.addReaction(sUptakeHPHB, -(1 - YSP));
		rsOxygen.addReaction(growthHPHB, -(YSP / YSX - 1));
		rsOxygen.addReaction(phbConsumption, -(1 - YPX));
		rsOxygen.addReaction(sUptakeH2, -(1 - YSP));
		rsOxygen.addReaction(growthHEPS, -(YSP / YSX - 1));
		oxygen.setProcesses(rsOxygen);
		// ---- Add the species to system ----
		addBiomassSpecies(speciesHPHB);
		addBiomassSpecies(speciesHEPS);
		addSoluteSpecies(oxygen);
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
		int[] nCells = { initialCellNumber, initialCellNumber };
		inoculateRandomly(nCells);
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
		ApplicationComponent app = new CaseStudyPHBvsEPS();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// uncomment the following line for plot of produced biomass
		// app = new SeriesVizualizer(app, prod);
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// uncomment the following line for plot of total biomass in biofilm
		// app = new SeriesVizualizer(app, biomass);
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
			// initialize
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
			// wait for user to press start iteration
			// app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}