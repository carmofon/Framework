package nl.tudelft.bt.model.tests.eps;
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
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
/**
 * Simple 2D monospecies system, can be used to study the effect of diffusion
 * limitation on biofilm structure (change growth rate and density, for
 * instance)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DRG2SpeciesEPSTest extends ModelHandler {
	//output directory name
	private final static String outputDirectory = "D:\\results"
			+ "\\multigrid_boundary_test\\test\\";
	// geometry (2D or 3D)
	private static int geometry = 2;
	//
	// Solute species
	// O2
	private static final float diffusivityO2 = 8333333f; //[um2/h]
	private static float bulkConcentrationO2 = 5e-3f; //[gO2/l]
	// S
	private static final float diffusivityS = 4166666f; //[um2/h]
	private static final float KS = 4e-3f; //[gCOD/l]
	private static float inputConcentrationS = 3.09E-04f; //[gCOD/l]
	// N
	private static final float diffusivityN = 7083333.333f; //[um2/h]
	private static final float KN = 0.0015f; //[gN/l]
	private static float inputConcentrationN = 0.002112985f; //[gN/l]
	//
	// Particulate species
	//Heterotroph
	private static float uMaxH = 0.2499f; // [1/h]
	private static float densityH = 10.0f; // [gX/l]
	private static final float YS = 0.63f; // [gX/gCOD]
	private static final float KO2_H = 2e-4f; //[gO2/l]
	//Nitrifier
	private static float uMaxN = 0.05775f; // [1/h]
	private static float densityN = densityH; // [gX/l]
	private static final float YN = 0.063f; // [gX/gN]
	private static final float KO2_N = 0.0005f; //[gO2/l]
	//EPS
	private static float densityEPS = densityH / 6f; // [gX/l]
	private static final float YEPS_H = 1 / 1.4f; // [gEPS/gH]
	//
	//Computational parameters
	private static float systemSize = 2000; // [um]
	private static float relativeMaximumRadius = 0.008f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static final float minimumMassRatio = 0.001f;
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	// reactor operation parameters
	private static float relativeBoundaryLayer = 0f;
	private static float maximumThickness = 540; // [um]
	private static int initialCellNumber = 16; // particles of each specie
	private static float reactorVolume = 1.25e15f; //[um^3]
	private static float carrierArea = 1E+11f; //[um^2]
	private static float residenceTime = 3f; //[h]
	private static float rdetach = 5e-5f;
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the solutes
		// oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", diffusivityO2);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				bulkConcentrationO2));
		// carbon source
		SoluteSpecies carbonSource = new SoluteSpecies("carbonSource",
				diffusivityS);
		carbonSource.setBulkConcentration(new ConstantBulkConcentration(
				inputConcentrationS));
		// nitrogen source
		SoluteSpecies nitrogenSource = new SoluteSpecies("nitrogenSource",
				diffusivityN);
		nitrogenSource.setBulkConcentration(new ConstantBulkConcentration(
				inputConcentrationN));
		//
		// create the species
		//heterotrophs
		ParticulateSpecies activeMassH = new ParticulateSpecies(
				"heterotrophActiveMass", densityH, Color.gray);
		ParticulateSpecies eps = new ParticulateSpecies("EPS", densityEPS, Color.gray);
		// array of fixed species that constitute the heterotroph
		ParticulateSpecies[] fH = {activeMassH, eps};
		float[] fractionalVolumeCompositionH = {0.5f, 0.5f};
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", fH,
				fractionalVolumeCompositionH);
		heterotroph.setActiveMass(activeMassH);
		heterotroph.setEpsMass(eps);
		//heterotrophs
		ParticulateSpecies activeMassN = new ParticulateSpecies(
				"nitrifierActiveMass", densityN, Color.gray);
		// array of fixed species that constitute the nitrifier
		ParticulateSpecies[] fN = {activeMassN};
		float[] fractionalVolumeCompositionN = {1.0f};
		// create the biomass species
		BiomassSpecies nitrifier = new BiomassSpecies("nitrifier", fN,
				fractionalVolumeCompositionN);
		//
		// create the reactions
		//heterotroph growth
		Reaction growthH = new Reaction("growthH",
				activeMassH, uMaxH, 1);
		growthH.addFactor(new MonodTwoSubstrates(oxygen, KO2_H, carbonSource,
				KS));
		//nitrifier growth
		Reaction growthN = new Reaction("growthN",
				activeMassN, uMaxN, 1);
		growthN.addFactor(new MonodTwoSubstrates(oxygen, KO2_N, nitrogenSource,
				KN));
		//
		// assign reaction to the species through a ReactionStoichiometry
		//heterotroph acive mass
		NetReaction rsXH = new NetReaction(1);
		rsXH.addReaction(growthH, 1.0f);
		activeMassH.setProcesses(rsXH);
		//eps
		NetReaction rsEPS = new NetReaction(1);
		rsEPS.addReaction(growthH, YEPS_H);
		eps.setProcesses(rsEPS);
		//nitrifier acive mass
		NetReaction rsXN = new NetReaction(1);
		rsXN.addReaction(growthN, 1.0f);
		activeMassN.setProcesses(rsXN);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsO2 = new NetReaction(2);
		rsO2.addReaction(growthH, -(1 - YS) / YS);
		rsO2.addReaction(growthN, -(4.57f - YN) / YN);
		oxygen.setProcesses(rsO2);
		// carbon source
		NetReaction rsS = new NetReaction(1);
		rsS.addReaction(growthH, -1 / YS);
		carbonSource.setProcesses(rsS);
		// nitrogen source
		NetReaction rsN = new NetReaction(1);
		rsN.addReaction(growthN, -1 / YN);
		nitrogenSource.setProcesses(rsN);
		// add the species to system
		addBiomassSpecies(heterotroph);
		addBiomassSpecies(nitrifier);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(carbonSource);
		addSoluteSpecies(nitrogenSource);
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
			int[] nCells = {initialCellNumber, initialCellNumber};
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
		// To use a detachment function set here
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
		//		try {
		//			Model.model().setMaximumBiofilmHeight(maximumThickness);
		//		} catch (InvalidValueException e) {
		//			throw new ModelRuntimeException(e);
		//		}
	}
	/**
	 * Simulation, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(10, 100);
		ApplicationComponent app = new DRG2SpeciesEPSTest();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		app = new DetachedBiomassVizualizer(app);
		// bulk concentrations
		app = new BulkConcentrationVizualizer(app);
		// finally, the controller must be the last decorator to add
		app = new VizualModelControler(app);
		float dX = systemSize * (1 - relativeBoundaryLayer);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			// set the global mass balance parameters
			app.setReactorParameters(residenceTime, carrierArea, reactorVolume);
			//initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			//app.addStateWritter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			//app.addStateWritter(new SimulationParametersWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer().getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			// add bulk concentrations of all solutes as variable series
			//
			app.addStateWriter(spw);
			//
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// add solute global rates to write
			for (Iterator iter = Model.model().getSoluteSpecies().iterator(); iter
					.hasNext();) {
				SoluteSpecies s = (SoluteSpecies) iter.next();
				spw.addSeries(s.getBulkConcentrationSeries());
			}
			// add particulate global masses to write
			for (Iterator iter = Model.model().getParticulateSpecies()
					.iterator(); iter.hasNext();) {
				ParticulateSpecies s = (ParticulateSpecies) iter.next();
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
			app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		//TODO remove this line
		System.out.println("done.");
	}
}