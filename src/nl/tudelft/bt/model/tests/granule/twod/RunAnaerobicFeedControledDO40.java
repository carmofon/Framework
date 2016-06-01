package nl.tudelft.bt.model.tests.granule.twod;

import java.awt.Color;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.work.granule.ItermittentBulkConcentrationControlled;

/**
 * Simulates the growth of a granule in 2D with fixed radius Note: um
 * representes micro-meter
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class RunAnaerobicFeedControledDO40 extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	//output directory name
	protected static String outputDirectory = "/Users/jxavier/results/lixoControled";
	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	//Reaction parameters
	//ammonium Oxydizers (XNH)
	protected static float uMaxXNH = 0.0854f; //[gCOD-XNH/gCOD-XNH/h]
	private static float K_XNH_O2 = 0.6e-3f; //[gO2/L]
	private static float K_XNH_NH4 = 2.4e-3f; //[gN/L]
	private static float Y_XNH_NH4 = 0.150f; // [gCOD-XNH/gN]
	//Nitrite oxydizers (XNO)
	protected static float uMaxNO = 0.0604f; //[gCOD-XNO/gCOD-XNO/h]
	private static float K_XNO_O2 = 2.2e-3f; //[gO2/L]
	private static float K_XNO_NO2 = 5.5e-3f; //[gN/L]
	private static float Y_XNO_NO2 = 0.041f; // [gCOD-XNO/gN]
	//Heterotrophs (XH)
	protected static float qSMaxXH = 0.952f; //[gCOD-XH/gCOD-XH/h]
	private static float K_XH_O2 = 3.5e-4f; //[gO2/l]
	private static float K_XH_S = 4e-3f; //[gCOD/l]
	private static float Y_H_S = 0.495f; //[gCOD-XH/gCOD-S]
	private static float fPHBmaxXH = 0.33f; //[gCOD-PHB/gCOD-XH]
	private static float Y_PHB_S_XH = 0.667f; //[gCOD-PHB/gCOD-S]
	protected static float kPhbXH = 0.15f; //[gCOH-PHB/gCOD-PHB/h]
	private static float Y_XH_PHB = 0.668f; //[gCOD-H/gCOD-PHB]
	//Phospate accumulating organisms (XPAO)
	protected static float qSMaxPAO = 0.4029f; //[gCOD-S/gCOD-XPAO/h]
	private static float KPAOO2 = 7e-4f; //[gO2/L]
	private static float KPAOS = 4e-3f; //[gCOD/L]
	private static float KPAOPP = 0.01e-3f; //[gP/L]
	private static float KPAOGLY = 0.01e-3f; //[gCOD-Gly/L]
	private static float Y_PHB_S_of_PAO = 1.5f; //[gCOD-PHB/gCOD-S]
	private static float Y_PO4_S = 0.3f; //[gP/gCOD-S]
	private static float fPHBmaxPAO = 1.0f; //[gCOD-PHB/gCOD-XPAO]
	//protected static float kPhbPAO = 0.3f; //[gCOH-PHB/gCOD-PHB/h]
	protected static float kPhbPAO = 0.15f; //[gCOH-PHB/gCOD-PHB/h]
	private static float YPhbAerobicPAO = 1.39f; //[gCOD-PHB/gCOD-XPAO]
	//private static float YPhbAerobicPAO = 1.5f; //[gCOD-PHB/gCOD-XPAO]
	private static float YPhbAnoxicPAO = 1.7f; //[gCOD-PHB/gCOD-XPAO]
	private static float KPAOPO4 = 3.1e-3f; //[gP/L]
	protected static float kPPPAO = 0.0187f; //[gP/gCOD-PHB/h]
	private static float YPpAerobicPAO = 4.37f; //[gP/gCOD-PAO]
	private static float YPpAnoxicPAO = 2.92f; //[gP/gCOD]
	private static float fPPmaxPAO = 0.5f; //[gP/gCOD-XPAO]
	private static float Y_GLY_PHB = 0.73f; //[gCOD-Gly/gCOD-PHB]
	protected static float kGlyPAO = 0.0454f; //[gCOD-Gly/gCOD-PHB/h]
	private static float YGlyAerobicPAO = 1.12f; //[gCOD-Gly/gCOD-XPAO]
	private static float YGlyAnoxicPAO = 1.15f; //[gCOD-Gly/gCOD]
	private static float fGLYmaxPAO = 0.45f; //[gCOD-Gly/gCOD-PAO]
	//Common to all species
	protected static float kDecay = 0.0033f; //[gCOD/gCOD/h]
	private static float fIDecay = 0.4f; //[gCOD-I/gCOD-X]
	private static float iNBM = 0.07f; //[gN/gCOD-X]
	private static float iNXI = 0.02f; //[gN/gCOD-I]
	private static float KNH4 = 0.05e-3f; //[gN/L]
	private static float KNO3 = 1e-3f; //[gN/L]
	private static float eta = 0.2f;

	// geometry (default is 3D) - change value to 2 for 2D
	protected static int geometry = 2;

	//Reactor properties
	protected static float reactorVolume = 3; // [L]
	protected static float feedFraction = 0.5f; // [dimensionless]

	//duration of the SBR cycle
	protected static float cycleTime = 3; //[h]
	//The SBR sycle
	private static SbrBulkConcentration.SbrCycle _cycle = new SbrBulkConcentration.SbrCycle(
			cycleTime);
	private static int _nCycles = 400; // write full data concerning cycle every
	// _nCycles

	//Feed composition
	//protected static float oxygenBulkConcentration = 10e-3f; //[gO2/L] DO 100
	protected static float oxygenBulkConcentration = 4e-3f; //[gO2/L] DO 20

	protected static float ammoniumFeedConcentration = 71e-3f; //[gN/L]

	protected static float nitriteFeedConcentration = 0; //[gN/L]

	protected static float nitrateFeedConcentration = 0; //[gN/L]

	protected static float substrateFeedConcentration = 396e-3f; //[gCOD/L]

	protected static float phosphateBulkConcentration = 20e-3f; //[gP/L]

	private static float oxygenDiffusivity = 8.3e6f; //[um^2/h]

	private static float substrateDiffusivity = 4e6f; //[um2/h]

	private static float phosphateDiffusivity = 4e6f; //[um2/h] ASSUMED!

	private static float ammoniumDiffusivity = 7.083e6f; //[um^2/h]

	private static float nitrateDiffusivity = 6.6667e6f; //[um^2/h]

	private static float nitriteDiffusivity = 6.6667e6f; //[um^2/h]

	//Leave these values
	protected static float nConcentrationsPrecision = 0.3e-3f; //[gN/L]

	protected static float sConcentrationsPrecision = 4e-3f; //[gCOD/L]

	protected static float pConcentrationsPrecision = 0.15e-3f; //[gP/L]

	//
	//Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; //[gCOD-H/L]

	protected static float specificMassPolymers = specificMassBiomass * 1e5f; //[gCOD-PHB/L]

	// Computation parameters
	//Size of computational volume (size of size of square)
	protected static float systemSize = 1600; // [um]

	//relativeMaximumRadius defines the maximum radius of the biomass particles
	//in relation to the system size
	//the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.007f;

	//Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 10 / systemSize;

	protected static float maximumGranuleRadius = 550; //[um]

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side
	//Don't change this

	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	//Don't change this

	//detachment rate
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]
	//leave at zero to form round granules

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumberXNH = 20;
	protected static int initialParticleNumberXNO = 20;
	protected static int initialParticleNumberXH = 20;
	protected static int initialParticleNumberXPAO = 20;

	//iteration finish time
	protected static float simulationFinishTime = 8000f; //[h]

	//outpute (write results to file) every:
	protected static float outputEvery = 15.0f; //[h]

	//Computational volume multiplier
	protected static float nComp = 25.93f * 6.7e5f; //[dimensionless]
	//factor 25.93f makes substrate loading per biomass the same as in real
	// reactor

	///END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		//1. Create the solutes
		//ammonium (NH3)
		SoluteSpecies ammonium = new SoluteSpecies("ammonium",
				ammoniumDiffusivity);
		//set up the bulk concentration
		ammonium.setBulkConcentration(new SbrBulkConcentration(
				ammoniumFeedConcentration, feedFraction,
				nConcentrationsPrecision, _cycle));
		//nitrite (NO2)
		SoluteSpecies nitrite = new SoluteSpecies("nitrite", nitriteDiffusivity);
		//set up the bulk concentration
		nitrite.setBulkConcentration(new SbrBulkConcentration(
				nitriteFeedConcentration, feedFraction,
				nConcentrationsPrecision, _cycle));
		//nitrite (NO3)
		SoluteSpecies nitrate = new SoluteSpecies("nitrate", nitrateDiffusivity);
		//set up the bulk concentration
		nitrate.setBulkConcentration(new SbrBulkConcentration(
				nitriteFeedConcentration, feedFraction,
				nConcentrationsPrecision, _cycle));
		//substrate (S)
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		//oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		//set up the bulk concentration for oxygen, uncomment one of these 3:
		//full aeration
		//oxygen.setBulkConcentration(new ConstantBulkConcentration(
		//		oxygenBulkConcentration));
		//Anaerobic feed phase
		//oxygen.setBulkConcentration(new ItermittentBulkConcentration(
		//		oxygenBulkConcentration, 0, 2f, 1f));
		//Anaerobic feed phase with controlled shutdown
		oxygen.setBulkConcentration(new ItermittentBulkConcentrationControlled(
				oxygenBulkConcentration, 0, 2f, 1f, ammonium, nitrite,
				nConcentrationsPrecision));
		//set up the bulk concentration
		substrate.setBulkConcentration(new SbrBulkConcentration(
				substrateFeedConcentration, feedFraction,
				sConcentrationsPrecision, _cycle));
		//phosphate (PO4)
		SoluteSpecies phosphate = new SoluteSpecies("phosphate",
				phosphateDiffusivity);
		//set up the bulk concentration
		phosphate.setBulkConcentration(new SbrBulkConcentration(
				phosphateBulkConcentration, feedFraction,
				pConcentrationsPrecision, _cycle));
		//2. Create the particulate species (solids)
		//NH
		ParticulateSpecies activeNH = new ParticulateSpecies("activeNH",
				specificMassBiomass, Color.blue);
		ParticulateSpecies inertNH = new ParticulateSpecies("inertNH",
				specificMassBiomass, Color.blue);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spNH = {activeNH, inertNH};
		float[] fractionalVolumeCompositionNH = {1.0f, 0};
		//3. Create the biomass species
		BiomassSpecies speciesNH = new BiomassSpecies("speciesNH", spNH,
				fractionalVolumeCompositionNH);
		speciesNH.setActiveMass(activeNH);
		speciesNH.setInertMass(inertNH);
		//NHO
		ParticulateSpecies activeNO = new ParticulateSpecies("activeNO",
				specificMassBiomass, Color.red);
		ParticulateSpecies inertNO = new ParticulateSpecies("inertNO",
				specificMassBiomass, Color.red);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spNO = {activeNO, inertNO};
		float[] fractionalVolumeCompositionNO = {1.0f, 0};
		//3. Create the biomass species
		BiomassSpecies speciesNO = new BiomassSpecies("speciesNO", spNO,
				fractionalVolumeCompositionNO);
		speciesNO.setActiveMass(activeNO);
		speciesNO.setInertMass(inertNO);
		//Heterotroph (H)
		ParticulateSpecies activeH = new ParticulateSpecies("activeH",
				specificMassBiomass, Color.orange);
		ParticulateSpecies inertH = new ParticulateSpecies("inertH",
				specificMassBiomass, Color.orange);
		ParticulateSpecies phbH = new ParticulateSpecies("phbH",
				specificMassPolymers, Color.orange);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spH = {activeH, phbH, inertH};
		float[] fractionalVolumeCompositionH = {1, 0, 0};
		//3. Create the biomass species
		BiomassSpecies speciesH = new BiomassSpecies("speciesH", spH,
				fractionalVolumeCompositionH);
		speciesH.setActiveMass(activeH);
		speciesH.setInertMass(inertH);
		//Pao (PAO)
		Color paoColor = Color.cyan;
		ParticulateSpecies activePAO = new ParticulateSpecies("activePAO",
				specificMassBiomass, paoColor);
		ParticulateSpecies inertPAO = new ParticulateSpecies("inertPAO",
				specificMassBiomass, paoColor);
		ParticulateSpecies phbPAO = new ParticulateSpecies("phbPAO",
				specificMassPolymers, paoColor);
		ParticulateSpecies polypPAO = new ParticulateSpecies("polypPAO",
				specificMassPolymers, paoColor);
		ParticulateSpecies glycogenPAO = new ParticulateSpecies("glycogenPAO",
				specificMassPolymers, paoColor);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spPAO = {activePAO, phbPAO, polypPAO, glycogenPAO,
				inertPAO};
		float initialFractionPAO = 0.4f / specificMassBiomass;
		float initialFractionPhb = 0.4f / specificMassPolymers;
		float initialFractionPolyP = 0.1f / specificMassPolymers;
		float initialFractionGlyP = 0.1f / specificMassPolymers;
		float total = initialFractionPAO + initialFractionPolyP
				+ initialFractionGlyP + initialFractionPhb;
		float[] fractionalVolumeCompositionPAO = {initialFractionPAO / total,
				initialFractionPhb / total, initialFractionPolyP / total,
				initialFractionGlyP / total, 0};
		//3. Create the biomass species
		BiomassSpecies speciesPAO = new BiomassSpecies("speciesPAO", spPAO,
				fractionalVolumeCompositionPAO);
		speciesPAO.setActiveMass(activePAO);
		speciesPAO.setInertMass(inertPAO);
		//4. Create the Reaction factors, Monod and inhibition coefficients
		//		ProcessFactor mS = new Saturation(oxygen, KO);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//for NH
		ProcessFactor mNHO2 = new Saturation(oxygen, K_XNH_O2);
		ProcessFactor iNHO2 = new Inhibition(oxygen, K_XNH_O2);
		ProcessFactor mNHNH4 = new Saturation(ammonium, K_XNH_NH4);
		//for NO
		ProcessFactor mNOO2 = new Saturation(oxygen, K_XNO_O2);
		ProcessFactor mNONO2 = new Saturation(nitrite, K_XNO_NO2);
		ProcessFactor iNOO2 = new Inhibition(oxygen, K_XNO_O2);
		//for H
		ProcessFactor mHO2 = new Saturation(oxygen, K_XH_O2);
		ProcessFactor mHS = new Saturation(substrate, K_XH_S);
		ProcessFactor iHS = new Inhibition(substrate, K_XH_S);
		ProcessFactor iHO2 = new Inhibition(oxygen, K_XH_O2);
		ProcessFactor iPHB = new InhibitionFromFractionCapacity(phbH, activeH,
				0.01f, fPHBmaxXH);
		//for PAO
		ProcessFactor mPAOS = new Saturation(substrate, KPAOS);
		ProcessFactor iPAOO2 = new Inhibition(oxygen, KPAOO2);
		ProcessFactor mPAOGLY = new Saturation(glycogenPAO, KPAOGLY);
		ProcessFactor mPAOPP = new Saturation(polypPAO, KPAOPP);
		ProcessFactor mPAOO2 = new Saturation(oxygen, KPAOO2);
		ProcessFactor mPAOPO4 = new Saturation(phosphate, KPAOPO4);
		ProcessFactor mf_PAOPHB = new SaturationFromFraction(phbPAO, activePAO,
				0.33f);
		ProcessFactor if_PAOPHB = new InhibitionFromFractionCapacity(phbPAO,
				activePAO, 0.01f, fPHBmaxPAO);
		ProcessFactor iPAOPP = new InhibitionFromFractionCapacity(polypPAO,
				activePAO, 0.01f, fPPmaxPAO);
		ProcessFactor iPAOGly = new InhibitionFromFractionCapacity(glycogenPAO,
				activePAO, 0.01f, fGLYmaxPAO);
		//for all organisms
		ProcessFactor mNO3 = new Saturation(nitrate, KNO3);
		ProcessFactor mNH4 = new Saturation(ammonium, KNH4);
		//5. Create the reactions
		//NH
		//growth NH
		Reaction growthNH = new Reaction("growthNH", activeNH, uMaxXNH, 2);
		growthNH.addFactor(mNHO2);
		growthNH.addFactor(mNHNH4);
		//NH decay
		Reaction decayNH = new Reaction("decayNH", activeNH, kDecay, 0);
		//NO
		//growth NO
		Reaction growthNO = new Reaction("growthNO", activeNO, uMaxNO, 2);
		growthNO.addFactor(mNOO2);
		growthNO.addFactor(mNONO2);
		//growthNO.addFactor(mNH4); //
		//NO decay
		Reaction decayNO = new Reaction("decayNO", activeNO, kDecay, 0);
		//H
		//aerobic substrate uptake H
		Reaction aerobicSUptakeH = new Reaction("aerobicSUptakeH", activeH,
				qSMaxXH, 3);
		aerobicSUptakeH.addFactor(mHO2);
		aerobicSUptakeH.addFactor(mHS);
		aerobicSUptakeH.addFactor(mNH4);
		//aerobic phb storage H
		Reaction aerobicPhbStorageH = new Reaction("aerobicPhbStorageH",
				activeH, qSMaxXH * Y_PHB_S_XH, 3);
		aerobicPhbStorageH.addFactor(mHO2);
		aerobicPhbStorageH.addFactor(iPHB);
		aerobicPhbStorageH.addFactor(mHS);
		//anoxic substrate uptake H
		Reaction anoxicSUptakeH = new Reaction("anoxicSUptakeH", activeH,
				qSMaxXH * eta, 4);
		anoxicSUptakeH.addFactor(iHO2);
		anoxicSUptakeH.addFactor(mNO3);
		anoxicSUptakeH.addFactor(mHS);
		anoxicSUptakeH.addFactor(mNH4);
		//anoxic phb storage H
		Reaction anoxicPhbStorageH = new Reaction("anoxicPhbStorageH", activeH,
				qSMaxXH * Y_PHB_S_XH * eta, 4);
		anoxicPhbStorageH.addFactor(iHO2);
		anoxicPhbStorageH.addFactor(mNO3);
		anoxicPhbStorageH.addFactor(iPHB);
		anoxicPhbStorageH.addFactor(mHS);
		//Phbconsumption
		Reaction aerobicPhbConsumptionH = new Reaction(
				"aerobicPhbConsumptionH", phbH, kPhbXH, 3);
		aerobicPhbConsumptionH.addFactor(mHO2);
		aerobicPhbConsumptionH.addFactor(iHS);
		aerobicPhbConsumptionH.addFactor(mNH4);
		//aerobic substrate uptake H
		//anoxicPhbConsumptionH
		Reaction anoxicPhbConsumptionH = new Reaction("anoxicPhbConsumptionH",
				phbH, kPhbXH * eta, 4);
		anoxicPhbConsumptionH.addFactor(iHO2);
		anoxicPhbConsumptionH.addFactor(iHS);
		anoxicPhbConsumptionH.addFactor(mNO3);
		anoxicPhbConsumptionH.addFactor(mNH4);
		//H decay
		Reaction decayH = new Reaction("decayH", activeH, kDecay, 0);
		//PAO
		//substrate uptake
		Reaction sUptakePAO = new Reaction("sUptakePAO", activePAO, qSMaxPAO, 4);
		sUptakePAO.addFactor(mPAOS);
		sUptakePAO.addFactor(if_PAOPHB);
		sUptakePAO.addFactor(mPAOGLY);
		sUptakePAO.addFactor(mPAOPP);
		//aerobic consumption of PHB
		Reaction aerobicPhbConsumptionPAO = new Reaction(
				"aerobicPhbConsumptionPAO", activePAO, kPhbPAO, 3);
		aerobicPhbConsumptionPAO.addFactor(mPAOO2);
		aerobicPhbConsumptionPAO.addFactor(mf_PAOPHB);
		aerobicPhbConsumptionPAO.addFactor(mNH4);
		//aerobic storage of PolyP
		Reaction aerobicStoragePolyPPAO = new Reaction(
				"aerobicStoragePolyPPAO", activePAO, kPPPAO, 3);
		aerobicStoragePolyPPAO.addFactor(mPAOO2);
		aerobicStoragePolyPPAO.addFactor(mPAOPO4);
		aerobicStoragePolyPPAO.addFactor(iPAOPP);
		//aerobic storage of glycogen
		Reaction aerobicStorageGlyPAO = new Reaction("aerobicStorageGlyPAO",
				activePAO, kGlyPAO, 3);
		aerobicStorageGlyPAO.addFactor(mPAOO2);
		aerobicStorageGlyPAO.addFactor(iPAOGly);
		aerobicStorageGlyPAO.addFactor(mf_PAOPHB);
		//anoxic consumption of PHB
		Reaction anoxicPhbConsumptionPAO = new Reaction(
				"anoxicPhbConsumptionPAO", activePAO, kPhbPAO * eta, 4);
		anoxicPhbConsumptionPAO.addFactor(iPAOO2);
		anoxicPhbConsumptionPAO.addFactor(mNO3);
		anoxicPhbConsumptionPAO.addFactor(mf_PAOPHB);
		anoxicPhbConsumptionPAO.addFactor(mNH4);
		//anoxic storage of PolyP
		Reaction anoxicStoragePolyPPAO = new Reaction("anoxicStoragePolyPPAO",
				activePAO, kPPPAO * eta, 4);
		anoxicStoragePolyPPAO.addFactor(iPAOO2);
		anoxicStoragePolyPPAO.addFactor(mNO3);
		anoxicStoragePolyPPAO.addFactor(mPAOPO4);
		anoxicStoragePolyPPAO.addFactor(iPAOPP);
		//anoxic storage of glycogen
		Reaction anoxicStorageGlyPAO = new Reaction("anoxicStorageGlyPAO",
				activePAO, kGlyPAO * eta, 4);
		anoxicStorageGlyPAO.addFactor(iPAOO2);
		anoxicStorageGlyPAO.addFactor(mNO3);
		anoxicStorageGlyPAO.addFactor(iPAOGly);
		anoxicStorageGlyPAO.addFactor(mf_PAOPHB);
		//decay
		Reaction decayPAO = new Reaction("decayPAO", activePAO, kDecay, 0);
		//
		//6. Assign reaction to the species through ReactionStoichiometries
		//active mass NH
		NetReaction rsNHactive = new NetReaction(2);
		rsNHactive.addReaction(growthNH, 1);
		rsNHactive.addReaction(decayNH, -1);
		activeNH.setProcesses(rsNHactive);
		//inert NH
		NetReaction rsInertNH = new NetReaction(1);
		rsInertNH.addReaction(decayNH, fIDecay);
		inertNH.setProcesses(rsInertNH);
		//active NO
		NetReaction rsActiveNO = new NetReaction(2);
		rsActiveNO.addReaction(growthNO, 1);
		rsActiveNO.addReaction(decayNO, -1);
		activeNO.setProcesses(rsActiveNO);
		//inert NO
		NetReaction rsInertNO = new NetReaction(1);
		rsInertNO.addReaction(decayNO, fIDecay);
		inertNO.setProcesses(rsInertNO);
		//active mass H
		NetReaction rsHactive = new NetReaction(7);
		rsHactive.addReaction(aerobicSUptakeH, Y_H_S);
		rsHactive.addReaction(anoxicSUptakeH, Y_H_S);
		rsHactive.addReaction(aerobicPhbStorageH, -Y_H_S / Y_PHB_S_XH);
		rsHactive.addReaction(anoxicPhbStorageH, -Y_H_S / Y_PHB_S_XH);
		rsHactive.addReaction(aerobicPhbConsumptionH, Y_XH_PHB);
		rsHactive.addReaction(anoxicPhbConsumptionH, Y_XH_PHB);
		rsHactive.addReaction(decayH, -1);
		activeH.setProcesses(rsHactive);
		//PHB in H
		NetReaction rsPhbH = new NetReaction(4);
		rsPhbH.addReaction(aerobicPhbStorageH, 1);
		rsPhbH.addReaction(anoxicPhbStorageH, 1);
		rsPhbH.addReaction(aerobicPhbConsumptionH, -1);
		rsPhbH.addReaction(anoxicPhbConsumptionH, -1);
		phbH.setProcesses(rsPhbH);
		//inert H
		NetReaction rsInertH = new NetReaction(1);
		rsInertH.addReaction(decayH, fIDecay);
		inertH.setProcesses(rsInertH);
		// active mass PAO
		NetReaction rsPaoActive = new NetReaction(7);
		rsPaoActive.addReaction(aerobicPhbConsumptionPAO, 1 / YPhbAerobicPAO);
		rsPaoActive.addReaction(aerobicStorageGlyPAO, -1 / YGlyAerobicPAO);
		rsPaoActive.addReaction(anoxicPhbConsumptionPAO, 1 / YPhbAnoxicPAO);
		rsPaoActive.addReaction(anoxicStorageGlyPAO, -1 / YGlyAnoxicPAO);
		rsPaoActive.addReaction(decayPAO, -1);
		rsPaoActive.addReaction(aerobicStoragePolyPPAO, -1 / YPpAerobicPAO);
		rsPaoActive.addReaction(anoxicStoragePolyPPAO, -1 / YPpAnoxicPAO);
		activePAO.setProcesses(rsPaoActive);
		//PHB in PAO
		NetReaction rsPhbPAO = new NetReaction(3);
		rsPhbPAO.addReaction(sUptakePAO, Y_PHB_S_of_PAO);
		rsPhbPAO.addReaction(aerobicPhbConsumptionPAO, -1);
		rsPhbPAO.addReaction(anoxicPhbConsumptionPAO, -1);
		phbPAO.setProcesses(rsPhbPAO);
		//PP in PAO
		NetReaction rsPpPAO = new NetReaction(3);
		rsPpPAO.addReaction(sUptakePAO, -Y_PO4_S);
		rsPpPAO.addReaction(aerobicStoragePolyPPAO, 1);
		rsPpPAO.addReaction(anoxicStoragePolyPPAO, 1);
		polypPAO.setProcesses(rsPpPAO);
		//Glycogen in PAO
		NetReaction rsGlyPAO = new NetReaction(3);
		rsGlyPAO.addReaction(sUptakePAO, 1 - Y_PHB_S_of_PAO);
		rsGlyPAO.addReaction(aerobicStorageGlyPAO, 1);
		rsGlyPAO.addReaction(anoxicStorageGlyPAO, 1);
		glycogenPAO.setProcesses(rsGlyPAO);
		//inerts in PAO
		NetReaction rsInertPAO = new NetReaction(1);
		rsInertPAO.addReaction(decayPAO, fIDecay);
		inertPAO.setProcesses(rsInertPAO);
		//assign reaction stoichiometry to the solutes
		//oxygen
		NetReaction rsOxygen = new NetReaction(8);
		rsOxygen.addReaction(growthNH, 1 - 3.43f / Y_XNH_NH4);
		rsOxygen.addReaction(growthNO, 1 - 1.14f / Y_XNO_NO2);
		rsOxygen.addReaction(aerobicSUptakeH, -(1 - Y_H_S));
		rsOxygen.addReaction(aerobicPhbStorageH, -(Y_H_S / Y_PHB_S_XH - 1));
		rsOxygen.addReaction(aerobicPhbConsumptionH, -(1 - Y_XH_PHB));
		rsOxygen.addReaction(aerobicPhbConsumptionPAO, 1 / YPhbAerobicPAO - 1);
		rsOxygen.addReaction(aerobicStorageGlyPAO, 1 - 1 / YGlyAerobicPAO);
		rsOxygen.addReaction(aerobicStoragePolyPPAO, -(1 / YPpAerobicPAO));
		oxygen.setProcesses(rsOxygen);
		//ammonium
		NetReaction rsAmmonium = new NetReaction(18);
		rsAmmonium.addReaction(growthNH, -1 / Y_XNH_NH4 - iNBM);
		rsAmmonium.addReaction(decayNH, iNBM - fIDecay * iNXI);
		rsAmmonium.addReaction(growthNO, -iNBM);
		rsAmmonium.addReaction(decayNO, iNBM - fIDecay * iNXI);
		rsAmmonium.addReaction(aerobicSUptakeH, -iNBM * Y_H_S);
		rsAmmonium.addReaction(aerobicPhbStorageH, iNBM * Y_H_S / Y_PHB_S_XH);
		rsAmmonium.addReaction(aerobicPhbConsumptionH, -iNBM * Y_PHB_S_XH);
		rsAmmonium.addReaction(anoxicSUptakeH, -iNBM * Y_H_S);
		rsAmmonium.addReaction(anoxicPhbStorageH, iNBM * Y_H_S / Y_PHB_S_XH);
		rsAmmonium.addReaction(anoxicPhbConsumptionH, -iNBM * Y_PHB_S_XH);
		rsAmmonium.addReaction(decayH, iNBM - fIDecay * iNXI);
		rsAmmonium
				.addReaction(aerobicPhbConsumptionPAO, -iNBM / YPhbAerobicPAO);
		rsAmmonium.addReaction(aerobicStoragePolyPPAO, iNBM / YPpAerobicPAO);
		rsAmmonium.addReaction(aerobicStorageGlyPAO, iNBM / YGlyAerobicPAO);
		rsAmmonium.addReaction(anoxicPhbConsumptionPAO, -iNBM / YPhbAnoxicPAO);
		rsAmmonium.addReaction(anoxicStoragePolyPPAO, iNBM / YPpAnoxicPAO);
		rsAmmonium.addReaction(anoxicStorageGlyPAO, iNBM / YGlyAnoxicPAO);
		rsAmmonium.addReaction(decayPAO, iNBM - fIDecay * iNXI);
		ammonium.setProcesses(rsAmmonium);
		//nitrite (NO2)
		NetReaction rsNitrite = new NetReaction(2);
		rsNitrite.addReaction(growthNH, 1 / Y_XNH_NH4);
		rsNitrite.addReaction(growthNO, -1 / Y_XNO_NO2);
		nitrite.setProcesses(rsNitrite);
		//nitrate (NO3)
		NetReaction rsNitrate = new NetReaction(7);
		rsNitrate.addReaction(growthNO, 1 / Y_XNO_NO2);
		rsNitrate.addReaction(anoxicSUptakeH, -(1 - Y_H_S) / 2.86f);
		rsNitrate.addReaction(anoxicPhbStorageH,
				-(Y_H_S / Y_PHB_S_XH - 1) / 2.86f);
		rsNitrate.addReaction(anoxicPhbConsumptionH, -(1 - Y_XH_PHB) / 2.86f);
		rsNitrate.addReaction(anoxicPhbConsumptionPAO,
				(1 / YPhbAnoxicPAO - 1) / 2.86f);
		rsNitrate.addReaction(anoxicStoragePolyPPAO,
				-(1 / 2.86f / YPpAnoxicPAO));
		rsNitrate.addReaction(anoxicStorageGlyPAO,
				(1 - 1 / YPhbAnoxicPAO) / 2.86f);
		nitrate.setProcesses(rsNitrate);
		//substrate (S)
		NetReaction rsSubstrate = new NetReaction(7);
		rsSubstrate.addReaction(aerobicSUptakeH, -1);
		rsSubstrate.addReaction(anoxicSUptakeH, -1);
		rsSubstrate.addReaction(sUptakePAO, -1);
		rsSubstrate.addReaction(decayNH, 1 - fIDecay);
		rsSubstrate.addReaction(decayNO, 1 - fIDecay);
		rsSubstrate.addReaction(decayH, 1 - fIDecay);
		rsSubstrate.addReaction(decayPAO, 1 - fIDecay);
		substrate.setProcesses(rsSubstrate);
		//phosphate (PO4)
		NetReaction rsPhosphate = new NetReaction(3);
		rsPhosphate.addReaction(sUptakePAO, Y_PO4_S);
		rsPhosphate.addReaction(aerobicStoragePolyPPAO, -1);
		rsPhosphate.addReaction(anoxicStoragePolyPPAO, -1);
		phosphate.setProcesses(rsPhosphate);
		//
		//7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesNH);
		addBiomassSpecies(speciesNO);
		addBiomassSpecies(speciesH);
		addBiomassSpecies(speciesPAO);
		addSoluteSpecies(ammonium);
		addSoluteSpecies(nitrite);
		//NOTE: must be set after ammonium and nitrate!
		addSoluteSpecies(oxygen);
		addSoluteSpecies(nitrate);
		addSoluteSpecies(substrate);
		addSoluteSpecies(phosphate);
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
			int[] nCells = {initialParticleNumberXNH, initialParticleNumberXNO,
					initialParticleNumberXH, initialParticleNumberXPAO};
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
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(kdetach);
		setDetachmentHandler(df);
		// set the maximum granule radius
		try {
			Model.model().setMaximumBiofilmHeight(maximumGranuleRadius);
		} catch (InvalidValueException e) {
			System.out.println(e);
			System.exit(-1);
		}
	}

	/**
	 * Set the reactors dimensions
	 */
	private static void setTheReactorDimensions(ApplicationComponent app) {
		float rVIM = reactorVolume * 1e15f; //reactor volume in cubic
		// micrometer
		float carrierArea = nComp
				* Model.model().getComputationalVolumeCarrierArea();
		app.setReactorParameters(cycleTime, carrierArea, rVIM);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		//program is running localy (in my laptop)
		boolean locationRemote = false;
		if (args.length > 0) {
			outputDirectory = args[0];
			locationRemote = true;
		}
		//set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new RunAnaerobicFeedControledDO40();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VaribleSeries biovolume = new BiovolumeSeries();
		VaribleSeries[] runLengthSeries = {new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries()};
		// The following code will be omitted if no vizuals are desired
		if (!locationRemote) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// add vizualizer for solutes rates
			//app = new SoluteRateSeriesVizualizer(app);
			// the biovolume visualizer
			//		app = new SeriesVizualizer(app, biovolume);
			//runlength
			//app = new SeriesVizualizer(app, runLengthSeries, "Run Length");
			// detached biomass
			//		app = new DetachedBiomassVizualizer(app);
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
			// --- nothing to set in this case: constant bulk concentration
			//initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			//Pov witer is added twice
			app.addStateWriter(new SbrFullCycleStateWriterDecorator(
					new PovRayWriter(), _cycle, _nCycles));
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			app.addStateWriter(new SbrFullCycleStateWriterDecorator(
					new SoluteConcentrationWriter(), _cycle, _nCycles));
			app.addStateWriter(new SbrFullCycleStateWriterDecorator(
					new SolidsConcentrationWriter(), _cycle, _nCycles));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new SbrFullCycleStateWriterDecorator(
					new ParticlePositionWriter(), _cycle, _nCycles));
			//app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(runLengthSeries[0]);
			spw.addSeries(runLengthSeries[1]);
			spw.addSeries(runLengthSeries[2]);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// intialize the reactor dimensions
			setTheReactorDimensions(app);
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
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
			//app.waitForStartIteratingRequest();
			app.startIterating();
		} catch (Exception e1) {
			app.forceWriteState();
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
	}
}