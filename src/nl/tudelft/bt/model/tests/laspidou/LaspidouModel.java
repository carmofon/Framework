package nl.tudelft.bt.model.tests.laspidou;
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
 * Example for model description paper
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class LaspidouModel extends ModelHandler {
	// geometry (2D or 3D)
	private static int geometry = 2;
	//
	//Solute species
	//Oxygen
	private static float oxygenBulkConcentration = 9.2f; //[gO/m3]
	private static float oxygenDiffusivity = 0.000173f; //[m2/d]
	private static float KO = 0.1f; //[gO/m3]
	//Substrate
	private static float substrateBulkConcentration = 500f; //[g-CODS/m3]
	private static float substrateDiffusivity = 0.000138f; //[m2/d]
	private static float KS = 20f; //[g-CODS/m3]
	//UAP
	private static float uapBulkConcentration = 0f; //[g-CODS/m3]
	private static float uapDiffusivity = 0.000138f; //[m2/d]
	private static float Kuap = 100f; //[g-CODuap/m3]
	//BAP
	private static float bapBulkConcentration = 0f; //[g-CODS/m3]
	private static float bapDiffusivity = 0.000138f; //[m2/d]
	private static float Kbap = 85f; //[g-CODbap/m3]
	//Particulate species
	// heterotroph
	private static float densityH = 70000f; //[gCOD_X/m3]
	//EPS
	private static float densityEPS = 200000f; //[gCOD_EPS/m3]
	//Inert
	private static float densityI = 220000; //[gCOD_I/m3]
	//
	//Yield coefficients
	private static float YES = 0.18f; //[gCOD_EPS/gCOS_S]
	private static float YUS = 0.05f; //[gCOD_UAP/gCOS_S]
	private static float YXS = 0.34f; //[gCOD_X/gCOS_S]
	private static float YXP = 0.45f; //[gCOD_X/gCOS_P]
	private static float fd = 0.8f; //[decay fraction]
	//
	// Processes
	//Acetate uptake
	private static float qMaxS = 28.5f; //[gCOD_S/gCOD_X/d]
	//UAP uptake
	private static float qMaxUap = 1.8f; //[gCOD_UAP/gCOD_X/d]
	//BAP uptake
	private static float qMaxBap = 0.1f; //[gCOD_BAP/gCOD_X/d]
	//EPS hydrolysis
	private static float k_hyd = 0.17f; //[gCOD_EPS/gCOD_EPS/d]
	//decay
	private static float kDecay = 0.4f; //[gCOD_X/gCOD_X/h]
	//
	// Computation parameters
	private static float systemSize = 0.001f; // [m]
	private static float relativeMaximumRadius = 0.006f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static float relativeBoundaryLayer = 0.1f;
	// other model parameters
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	private static float rdetach = 1.4e8f; //detachment coef for computed prof
	//private static float rdetach = 5e9f; //for imposed O2 profiles
	// detachment constant[g/l/h]
	private static int initialCellNumber = 83;
	//private static int initialCellNumber = 31;
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
		//substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		//uap
		SoluteSpecies uap = new SoluteSpecies("uap", uapDiffusivity);
		uap.setBulkConcentration(new ConstantBulkConcentration(
				uapBulkConcentration));
		//bap
		SoluteSpecies bap = new SoluteSpecies("bap", bapDiffusivity);
		bap.setBulkConcentration(new ConstantBulkConcentration(
				bapBulkConcentration));
		// create the particulates
		//H active mass
		ParticulateSpecies activeH = new ParticulateSpecies(
				"heterotrophActiveMass", densityH, Color.blue);
		//EPS
		ParticulateSpecies eps = new ParticulateSpecies("eps", densityEPS,
				Color.yellow);
		//Inert
		ParticulateSpecies inert = new ParticulateSpecies("inert", densityI,
				Color.gray);
		// array of fixed species that constitute H2
		ParticulateSpecies[] spH = {activeH, eps, inert};
		float[] fractionalVolumeCompositionH1 = {1.0f, 0, 0};
		// create the biomass species
		BiomassSpecies speciesH = new BiomassSpecies("heterotroph", spH,
				fractionalVolumeCompositionH1);
		speciesH.setInertMass(inert);
		speciesH.setEpsMass(eps);
		//Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mO = new Saturation(oxygen, KO);
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor mUap = new Saturation(uap, Kuap);
		ProcessFactor mBap = new Saturation(bap, Kbap);
		// create the reactions
		//acetate uptake H1
		Reaction substrateUptake = new Reaction(
				"substrateUptake", activeH, qMaxS, 2);
		substrateUptake.addFactor(mO);
		substrateUptake.addFactor(mS);
		//uap uptake
		Reaction uapUptake = new Reaction("uapUptake",
				activeH, qMaxUap, 2);
		uapUptake.addFactor(mO);
		uapUptake.addFactor(mUap);
		//bap uptake
		Reaction bapUptake = new Reaction("bapUptake",
				activeH, qMaxBap, 2);
		bapUptake.addFactor(mO);
		bapUptake.addFactor(mBap);
		//eps hydrolysis
		Reaction epsHydrolysis = new Reaction(
				"epsHydrolysis", activeH, k_hyd, 0);
		//decay hydrolysis
		Reaction decay = new Reaction("decay", activeH,
				kDecay, 0);
		// assign reaction to the species through ReactionStoichiometries
		//H1 active mass
		NetReaction rsActiveH = new NetReaction(4);
		rsActiveH.addReaction(substrateUptake, YXS * (1 - YUS - YES));
		rsActiveH.addReaction(uapUptake, YXP);
		rsActiveH.addReaction(bapUptake, YXP);
		rsActiveH.addReaction(decay, -1);
		activeH.setProcesses(rsActiveH);
		//uap
		NetReaction rsUap = new NetReaction(2);
		rsUap.addReaction(substrateUptake, YUS);
		rsUap.addReaction(uapUptake, -1);
		uap.setProcesses(rsUap);
		//bap
		NetReaction rsBap = new NetReaction(2);
		rsBap.addReaction(epsHydrolysis, 1);
		rsBap.addReaction(bapUptake, -1);
		bap.setProcesses(rsBap);
		//eps
		NetReaction rsEps = new NetReaction(2);
		rsEps.addReaction(substrateUptake, YES);
		rsEps.addReaction(epsHydrolysis, -1);
		eps.setProcesses(rsEps);
		//inert
		NetReaction rsInert = new NetReaction(1);
		rsInert.addReaction(decay, (1 - fd));
		inert.setProcesses(rsInert);
		//substrate
		NetReaction rsS = new NetReaction(1);
		rsS.addReaction(substrateUptake, -1);
		substrate.setProcesses(rsS);
		//oxygen
		NetReaction rsO = new NetReaction(4);
		rsO.addReaction(substrateUptake, -(1 - YXS));
		rsO.addReaction(uapUptake, -(1 - YXP));
		rsO.addReaction(bapUptake, -(1 - YXP));
		rsO.addReaction(decay, -fd);
		oxygen.setProcesses(rsO);
		// add the species to system
		addBiomassSpecies(speciesH);
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
		//output directory name
		//String outputDirectory = "D:/results/"
		//	+ "laspidouModel/2D/computedO2_"+ rdetach +"/";
		String outputDirectory = "/Users/jxavier/model/lixo/";
		//
		MultigridVariable.setSteps(10, 100);
		ApplicationComponent app = new LaspidouModel();
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