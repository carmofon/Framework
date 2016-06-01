package nl.tudelft.bt.model.tests.airlift;
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
public class AirliftVariableS extends ModelHandler {
	//output directory name
	private static String outputDirectory = "/Users/jxavier/model/results/"
			+ "airliftStructure/variableS/sim3/";
	// geometry (2D or 3D)
	private static int geometry = 2;
	//
	//Solute species
	//Carbon source
	private static float carbonSourceBulkConcentration = 0.1f; //[gCOD_S/l]
	private static float carbonSourceDiffusivity = 4e6f; //[um2/h]
	private static float KS = 4e-3f; //[gO/l]
	//Oxygen
	private static float oxygenBulkConcentration = 4e-3f; //[gO/l]
	private static float oxygenDiffusivity = 8e6f; //[um2/h]
	private static float KO = 2e-4f; //[gO/l]
	//
	//Particulate species
	//Non PHB producing heterotroph H2
	private static float densityH = 200f; //[gCOD_X/l]
	//
	//Yield coefficients
	private static float YSX = 0.63f; //[gCOD_X/gCOS_S]
	//
	// Processes
	//Carbon source uptake
	private static float qMax = 1.0f; //[gCOD_S/gCOD_X/h]
	// Computation parameters
	private static float systemSize = 1000; // [um]
	private static float relativeMaximumRadius = 0.008f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static float relativeBoundaryLayer = 0.1f;
	// other model parameters
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	private static float rdetach = 3e-3f;
	// detachment constant[g/l/h]
	private static int initialCellNumber = 62;
	// reactor parameters
	private static float residenceTime = 2.5f; //[h]
	private static float carrierArea = 1e11f; //[um2]
	private static float reactorVolume = 1.25E15f; //[um3]
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the solutes
		//carbon source
		SoluteSpecies carbonSource = new SoluteSpecies("carbonSource",
				carbonSourceDiffusivity);
		carbonSource.setBulkConcentration(new DynamicBulkConcentrationImplicit(
				carbonSourceBulkConcentration));
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
		//Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mHAc = new Saturation(carbonSource, KS);
		ProcessFactor mO = new Saturation(oxygen, KO);
		// create the reactions
		//acetate uptake H1
		Reaction acetateUptakeH1 = new Reaction(
				"acetateUptakeH1", activeH1, qMax, 2);
		acetateUptakeH1.addFactor(mHAc);
		acetateUptakeH1.addFactor(mO);
		// assign reaction to the species through ReactionStoichiometries
		//H1 active mass
		NetReaction rsH1active = new NetReaction(1);
		rsH1active.addReaction(acetateUptakeH1, YSX);
		activeH1.setProcesses(rsH1active);
		// assign reaction stoichiometry to the solutes
		//acetate
		NetReaction rsAcetate = new NetReaction(1);
		rsAcetate.addReaction(acetateUptakeH1, -1);
		carbonSource.setProcesses(rsAcetate);
		//oxygen
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(acetateUptakeH1, -(1 - YSX));
		oxygen.setProcesses(rsOxygen);
		// add the species to system
		addBiomassSpecies(speciesH1);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(carbonSource);
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
			int[] nCells = {initialCellNumber};
			inoculateRandomly(nCells);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
			return;
		}
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#initializeDetachmentFunction()
	 */
	public void initializeDetachmentFunction() {
		// TODO set dtachment function here
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}
	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length == 1)
			outputDirectory = args[0];
		MultigridVariable.setSteps(10, 100);
		ApplicationComponent app = new AirliftVariableS();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		//app = new SeriesVizualizer(app, prod);
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		//app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		app = new DetachedBiomassVizualizer(app);
		// bulk concentrations
		app = new BulkConcentrationVizualizer(app);
		// finally, the controller must be the last decorator to add
		app = new VizualModelControler(app);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			//initialize
			app.initializeSystemSpace();
			// set the global mass balance parameters
			app.setReactorParameters(residenceTime, carrierArea, reactorVolume);
			//
			app.intializeStateWriters(outputDirectory);
			app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			app.addStateWriter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer().getTotalDetachedBiomassSeries());
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
		//TODO remove this line
		System.out.println("done.");
	}
}