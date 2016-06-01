package nl.tudelft.bt.model.tests.biofilmcontrol;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2EpsFractionDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.tests.biofilmcontrol.particleparser.ParticlePositionParser;

/**
 * Control simulations for the induced detachment paper. Growth with EPS
 * production carried out without interference of a PDP
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class CaseStudy3 extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/jxavier/results/2DSims/"
			+ "Stoichiometric";

	//	protected static String particleFile = "/Users/jxavier/results/"
	//			+ "control_rd_0.0001_cS_0.001_seed_101/particles/iteration4178.txt";
	protected static String particleFile = "/Users/jxavier/results/"
			+ "2DSims/control_rd_1.0E-4_cS_0.0010_seed_103/"
			+ "particles/iteration00000088.txt";

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

	// //Microbial parameters
	// Maximum specific growth rate
	protected static float uMax = 0.5f; // [gX/gX/h]

	// Substrate saturation constant
	private static float KS = 4e-4f; // [gS/L]

	// Yield of EPS on produced biomass
	private static float YEPS = 1.4f; // [gEPS/gX]

	// Yield of consumed substrate per biomass produced
	private static float YS = 1.63f; // [gS/gX]

	// Yield of EPS_star per EPS decayed
	private static float YEPS_star = 1; // [gEPS_star/gEPS]

	// Yield of PDP consumed per EPS decayed
	private static float YPDP = 0.01f; // [gPDP/gEPS]

	// diffusivity of solutes
	private static float substrateDiffusivity = 4.2e6f; // [um^2/h]

	// Enzyme parameters
	private static float epsDecayRate = 100f; // [gEPS/gPDP/h]

	private static float pdpDiffusivity = 2.1e4f; // [um^2/h]

	private static float applicationTime = 1000f; // [h]

	private static float beginApplicationAt = 1e20f; // [h]

	protected static float pdpBulkConcentration = 1f; // [gPDP/L]

	// Concentration of solutes
	protected static float substrateBulkConcentration = 1e-3f; // [gO/L]

	// Specific mass of particulates
	protected static float specificMassX = 200f; // [gX/L]

	protected static float specificMassEPS = specificMassX / 6f; // [gEPS/L]

	protected static float specificMassEPS_star = specificMassEPS; // [gEPS/L]

	// monod constant for EPD decay
	private static float KEPS = specificMassEPS * 0.1f; // [gEPS/gPDP/h]

	// detachment
	protected static float rdetach = 1e-4f;
	protected static float gamma = 10.0f;


	// uncomment
	// inoculation
	protected static int initialParticleNumber = 3200;

	// Numeric parameters
	protected static float systemSize = 1000; // [um]

	// grid resolution
	protected static int gridSide = 33; // multigrid grid side

	// comment
	// protected static int initialParticleNumber = 20;
	// protected static float systemSize = 500; // [um]
	// protected static int gridSide = 17; // multigrid grid side
	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 10 / systemSize;

	// protected static float relativeMaximumRadius = 14 / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = 100 / systemSize;

	// Shoving parameter
	protected static float kShov = 1.0f;

	// outpute (write results to file) every:
	protected static float outputEvery = 0.1f; // [h]

	// protected static float outputEvery = 0.5f; // [h]

	// eps and eps_star must be atributes since they are also used
	// todefine the detachment function
	private ParticulateSpecies _activeX;

	private ParticulateSpecies _eps;

	private ParticulateSpecies _epsStar;

	private BiomassSpecies _speciesX;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// //1. Create the solutes
		// substrate (S)
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		SoluteSpecies pdp = new SoluteSpecies("PDP", pdpDiffusivity);
		pdp.setBulkConcentration(new ConstantBulkConcentration(
				pdpBulkConcentration));
		// 2. Create the particulate species (soliids)
		// X active mass
		_activeX = new ParticulateSpecies("activeX", specificMassX, Color.blue);
		// EPS
		_eps = new ParticulateSpecies("EPS", specificMassEPS, Color.gray);
		// EPS*
		_epsStar = new ParticulateSpecies("EPS_star", specificMassEPS,
				Color.red);
		// array of fixed species that constitute speciesX (in this case
		// active mass, EPS and EPS*)
		ParticulateSpecies[] spX = {_activeX, _eps, _epsStar};
		// volunmetric fractional components
		float fEpsV = 1 / (1 + specificMassEPS / (YEPS * specificMassX));
		float[] fractionalVolumeCompositionH1 = {1 - fEpsV, fEpsV, 0};
		// 3. Create the biomass species
		_speciesX = new BiomassSpecies("speciesX", spX,
				fractionalVolumeCompositionH1);
		_speciesX.setActiveMass(_activeX);
		ParticulateSpecies[] eps_components = {_eps, _epsStar};
		_speciesX.setEpsMass(eps_components);
		// speciesX.getColorFromGrowth();
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor mEPS = new Saturation(_eps, KEPS);
		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", _activeX, uMax, 1);
		growth.addFactor(mS);
		// // EPS decay
		Reaction epsDecay = new Reaction("epsDecay", pdp, epsDecayRate, 1);
		epsDecay.addFactor(mEPS);
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(1);
		rsXactive.addReaction(growth, 1);
		_activeX.setProcesses(rsXactive);
		// EPS
		NetReaction rsEps = new NetReaction(2);
		rsEps.addReaction(growth, YEPS);
		rsEps.addReaction(epsDecay, -1);
		_eps.setProcesses(rsEps);
		// EPS*
		NetReaction rsEps_star = new NetReaction(1);
		rsEps_star.addReaction(epsDecay, YEPS_star);
		_epsStar.setProcesses(rsEps_star);
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsSubstrate = new NetReaction(1);
		rsSubstrate.addReaction(growth, -YS);
		substrate.setProcesses(rsSubstrate);
		// enzyme
		NetReaction rsEnzyme = new NetReaction(1);
		rsEnzyme.addReaction(epsDecay, -YPDP);
		pdp.setProcesses(rsEnzyme);
		// 7. add the solute species and the biomass
		// species (which contain the
		// particulate species) to system
		addBiomassSpecies(_speciesX);
		addSoluteSpecies(substrate);
		addSoluteSpecies(pdp);
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
		ParticlePositionParser parser = new ParticlePositionParser(
				particleFile, _speciesX, specificMassX, _activeX,
				specificMassEPS, _eps);
		parser.initializeSystemState(this);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		// DetachmentSpeedFunction df = new Height2EpsDetachment(rdetach,
		// specificMassEPS);
		// DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		DetachmentSpeedFunction df = new Height2EpsFractionDetachment(rdetach,
				_eps, _epsStar, gamma);
		setDetachmentHandler(df);
	}
	/*
	 * Overrides the method to check if the model should continue computation
	 */
	private final ArrayList _particleNumbers = new ArrayList();
	private final ArrayList _timeSeries = new ArrayList();
	public void computeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		super.computeSoluteConcentrations();
		//get the time
		float time = Model.model().getTime();

		// get the number of biomass partilces still in the system
		int nParticles = Model.model().getNumberOfParticles(_speciesX);
		_particleNumbers.add(new Integer(nParticles));
		_timeSeries.add(new Float(time));
		if (nParticles == 0) {
			//finish simulations if number of particles is 0
			Model.model().setFinishSimulation();
			System.out.println("No particles in system.");
		} else if (time >= 1) {
			// get the number of particles in previous iteration
			int nParticlesPrevious = ((Integer) _particleNumbers
					.get(_particleNumbers.size() - 2)).intValue();
			//measure the rate in particle change
			float rate = (nParticles - nParticlesPrevious)
					/ Model.model().getTimeStep();
			// determine a new timestep based on that rate
			float timeStep = 1;
			if (rate < 0) {
				timeStep = -0.1f * nParticles / rate;
				timeStep = (timeStep < 0.1f ? 0.1f : timeStep);
			}
			else if (rate > 0) {
				timeStep = 0.1f * nParticles / rate;
				timeStep = (timeStep < 1.0f ? 1.0f : timeStep);
			}
			// don't allow steps greater than 1 until 10 h
			if (time < 10.0f)
				timeStep = (timeStep > 1.0f ? 1.0f : timeStep);
			else 
				timeStep = (timeStep > 10.0f ? 10.0f : timeStep);				
			// set the time step
			Model.model().overrideTimeStep(timeStep);
			//if time is over 30, start checks:
			// register time every 10 h
			if ((time > 40.0f) & (time % 10.0f < 1.0f)) {
				//find the index of the latest point that it at least 10 h
				//before present time
				int indice10 = 0;
				//iterate over lists
				for (int i = 0; i < _timeSeries.size(); i++) {
					float newTime = ((Float) _timeSeries.get(i)).floatValue();
					if (time - newTime < 10.0f) {
						indice10 = i - 1;
						break;
					}
				}
				//if the number of particles has increased break
				int nParticlesNew = ((Integer) _particleNumbers.get(indice10))
						.intValue();
				if (nParticles > nParticlesNew) {
					// set iteration to finish
					Model.model().setFinishSimulation();
					System.out
							.println("long term fluctuation in particle number: "
									+ nParticles + " > " + nParticlesNew);
				}
			}
		}

	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my computer)
		// Arguments should be:
		// args[0] - parameter file name
		// args[1] - epsDecayRate * pdpBulkConcentration
		// args[2] - pdpDiffusivity / pdpDecayRate
		boolean locationRemote = false;
		if (args.length > 0) {
			locationRemote = true;
			geometry = 3;
			outputDirectory = "/home/jxavier/results/Stoichiometric";
			particleFile = args[0];
			epsDecayRate = Float.parseFloat(args[1]);
			YPDP = Float.parseFloat(args[2]);
			YEPS_star = Float.parseFloat(args[3]);
			System.out.println("Simulation started from file " + args[0]
					+ " with kdecEPS = " + args[1] + ", YPDP = " + args[2]
					+ ", YEPSstar = " + args[3]);
		}
		outputDirectory += "_kdecEPS_" + epsDecayRate + "_YPDP_" + YPDP
				+ "_YEPSstar_" + YEPS_star;
		// set multigrid variables:
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new CaseStudy3();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		if (!locationRemote) {
			// The following code will be omitted if no vizuals are desired
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
			// add visualizer for solutes rates
			app = new SoluteRateSeriesVizualizer(app);
			// add visualizer for bulk concentrations
			app = new BulkConcentrationVizualizer(app);
			// add the controler
			app = new VizualModelControler(app);
		}
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			// app.addTimedStateWriter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			try {
				// initialize
				app.initializeDiffusionReactionSystem(); // also innoculates
			} catch (ModelException e) {
				e.printStackTrace();
				System.exit(-1);
			}
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
			//
			Model.model().overrideTimeStep(outputEvery);
			// start iterating cycle
			//app.waitForStartIteratingRequest();
			app.startIterating();
		} catch (InterruptedException e1) {
			app.forceWriteState();
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
	}
}