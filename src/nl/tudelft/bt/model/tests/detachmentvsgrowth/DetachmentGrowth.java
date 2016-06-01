package nl.tudelft.bt.model.tests.detachmentvsgrowth;

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
 * Simulate the structure of biofilm in airlift reactor under detachment with
 * constant S concentration
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DetachmentGrowth extends ModelHandler {
	//output directory name
	protected static String outputDirectory = "~/results/lixo/";

	// geometry (2D or 3D)
	protected static int geometry = 2;

	//
	//Solute species
	//Oxygen
	protected static float oxygenBulkConcentration = 4e-3f; //[gO/l]

	private static float oxygenDiffusivity = 8e6f; //[um2/h]

	private static float KO = 3.5e-4f; //[gO/l]

	//
	//Particulate species
	//Non PHB producing heterotroph H2
	protected static float densityH = 200f; //[gCOD_X/l]

	//
	//Yield coefficients
	private static float YSX = 0.495f; //[gCOD_X/gCOS_S]

	//
	// Processes
	//Carbon source uptake
	protected static float qMax = 0.952f; //[gCOD_S/gCOD_X/h] Beun

	// Computation parameters
	protected static float systemSize = 2000; // [um]

	protected static float relativeMaximumRadius = 0.003f;

	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	protected static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 1.5e-5f;

	// detachment constant[g/l/h]
	protected static int initialCellNumber = 500;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the solutes
		//oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// create the particulates
		//H1 active mass
		ParticulateSpecies activeH1 = new ParticulateSpecies(
				"heterotroph1ActiveMass", densityH, Color.gray);
		// array of fixed species that constitute H1
		ParticulateSpecies[] spH1 = {activeH1};
		float[] fractionalVolumeCompositionH1 = {1.0f};
		// create the biomass species
		BiomassSpecies speciesH1 = new BiomassSpecies("heterotrophH1", spH1,
				fractionalVolumeCompositionH1);
		speciesH1.setActiveMass(activeH1);
		speciesH1.getColorFromGrowth();
		//Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO);
		// create the reactions
		//substrate uptake
		Reaction substrateUptake = new Reaction("acetateUptakeH1", activeH1,
				qMax, 1);
		substrateUptake.addFactor(mO);
		// assign reaction to the species through ReactionStoichiometries
		//active mass
		NetReaction rsH1active = new NetReaction(1);
		rsH1active.addReaction(substrateUptake, YSX);
		activeH1.setProcesses(rsH1active);
		// assign reaction stoichiometry to the solutes
		//oxygen
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(substrateUptake, -(1 - YSX));
		oxygen.setProcesses(rsOxygen);
		// add the species to system
		addBiomassSpecies(speciesH1);
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
		int[] nCells = {initialCellNumber};
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#initializeDetachmentFunction()
	 */
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param showVizuals
	 * @param args
	 */
	public void run(boolean showVizuals) {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = this;
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		if (showVizuals) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
			//app = new SeriesVizualizer(app, prod);
			//app = new SeriesVizualizer(app, biomass);
			// add vizualizer for solutes rates
			app = new SoluteRateSeriesVizualizer(app);
			// detached biomass
			app = new DetachedBiomassVizualizer(app);
			// bulk concentrations
			app = new BulkConcentrationVizualizer(app);
			// finally, the controller must be the last decorator to add
			app = new VizualModelControler(app);
		}
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
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			//app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
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
			app.writeState();
			app.startIterating();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		//TODO remove this line
		System.out.println("done.");
	}
}