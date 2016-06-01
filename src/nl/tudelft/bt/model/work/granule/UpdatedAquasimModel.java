package nl.tudelft.bt.model.work.granule;

import java.awt.Color;
import java.io.IOException;
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

/**
 * Simulates the growth of a granule in 2D with fixed radius. <br>
 * Note: um representes micro-meter
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class UpdatedAquasimModel extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory;

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// Reaction parameters
	// Ammonium oxydizers (XNH)
	protected static float uMaxXNH = 0.0167f; // [gCOD-XNH/gCOD-XNH/h] Merle

	private static float K_XNH_O2 = 0.3e-3f; // [gO2/L] Merle

	private static float K_XNH_NH4 = 2.4e-3f; // [gN/L]

	private static float Y_XNH_NH4 = 0.150f; // [gCOD-XNH/gN]

	// Nitrite oxydizers (XNO)
	protected static float uMaxNO = 0.0458f; // [gCOD-XNO/gCOD-XNO/h] Merle

	private static float K_XNO_O2 = 0.1e-3f; // [gO2/L] Merle

	private static float K_XNO_NO2 = 0.238e-3f; // [gN/L] Merle

	private static float Y_XNO_NO2 = 0.041f; // [gCOD-XNO/gN]

	// Heterotrophs: non-PAO (XH)
	protected static float uMaxH = 0.47f; // [gCOD-XH/gCOD-XH/h] xavier (j.

	// biofilms)

	private static float K_XH_O2 = 3.5e-4f; // [gO2/L] xavier (B&B)

	private static float K_XH_S = 4e-3f; // [gCOD-S/L] xavier (B&B)

	private static float Y_XH_S = 0.6645f; // [gCOD-XH/gS] xavier (j. biofilms)

	private static float K_XH_NO3 = 1e-3f; // [gN/L] same as for PAO

	// Phospate accumulating organisms (XPAO)
	protected static float qSMaxPAO = 0.40f; // [gCOD-S/gCOD-XPAO/h] Murnl

	private static float KPHA_P = 0.01f; // [gCOD-PHA/gCOD-XPAO]

	private static float KPP_P = 0.01f; // [gP/gCOD-XPAO]

	private static float KPAOO2 = 2e-4f; // [gO2/L]

	private static float KPAOS = 4e-3f; // [gCOD/L]


	private static float Y_PHB_S_of_PAO = 1.5f; // [gCOD-PHB/gCOD-S]

	// private static float Y_PO4_S = 0.3f; //[gP/gCOD-S]
	private static float Y_PO4_S = 0.5f; // [gP/gCOD-S] Merle

	private static float fPHBmaxPAO = 0.8f; // [gCOD-PHB/gCOD-XPAO] Merle

	protected static float kPhbPAO = 0.31f; // [gCOH-PHB/gCOD-PHB/h] Murnl

	private static float YPhbAerobicPAO = 1.39f; // [gCOD-PHB/gCOD-XPAO] Murnl

	private static float YPhbAnoxicNO3PAO = 1.7f; // [gCOD-PHB/gCOD-XPAO]

	private static float YPhbAnoxicNO2PAO = 1.7f; // [gCOD-PHB/gCOD-XPAO]

	private static float KPAOPO4 = 1e-6f; // [gP/L]

	private static float KPAOPO4_PP = 3.1e-3f; // [gP/L] Murnl

	private static float KfPHA_P = 0.33f; // [gCOD-PHA/gCOD-XPAO]

	protected static float kPPPAO = 0.019f; // [gP/gCOD-PAO/h] Murnl

	// private static float YPpAerobicPAO = 4.37f; //[gP/gCOD-XPAO]
	private static float YPpAerobicPAO = 4.42f; // [gP/gCOD-XPAO] Merle

	// private static float YPpAnoxicPAO = 2.92f; //[gP/gCOD]
	private static float YPpAnoxicNO3PAO = 3.02f; // [gP/gCOD] Merle

	private static float YPpAnoxicNO2PAO = 3.02f; // [gP/gCOD] Merle

	// private static float fPPmaxPAO = 0.5f; //[gP/gCOD-XPAO]
	private static float fPPmaxPAO = 0.65f; // [gP/gCOD-XPAO] Merle

	//protected static float kGlyPAO = 0.045f; // [gCOD-Gly/gCOD-PHB/h] Murnl
	protected static float kGlyPAO = 0.2f; // [gCOD-Gly/gCOD-PHB/h] TEST

	private static float YGlyAerobicPAO = 1.12f; // [gCOD-Gly/gCOD-XPAO] Murnleitner

	// private static float YGlyAnoxicPAO = 1.15f; //[gCOD-Gly/gCOD]
	private static float YGlyAnoxicNO3PAO = 1.18f; // [gCOD-Gly/gCOD] Merle

	private static float YGlyAnoxicNO2PAO = 1.18f; // [gCOD-Gly/gCOD] Merle

	// private static float fGLYmaxPAO = 0.45f; //[gCOD-Gly/gCOD-PAO]
	private static float fGLYmaxPAO = 0.5f; // [gCOD-Gly/gCOD-PAO] Merle

	// Common to all species
	protected static float kDecay = 0.0033f; // [gCOD/gCOD/h]

	// protected static float kDecay = 0.0167f; //[gCOD/gCOD/h] Merle

	// private static float fIDecay = 0.4f; //[gCOD-I/gCOD-X]
	private static float fIDecay = 1f; // [gCOD-I/gCOD-X] Merle

	private static float KPAONO3 = 1e-3f; // [gN/L]

	private static float KPAONO2 = 1e-3f; // [gN/L]

	private static float etaNO3 = 0.5f; // Mark

	private static float etaNO2 = 0.5f; // Mark

	// geometry (default is 3D) - change value to 2 for 2D
	protected static int geometry = 2;

	// Reactor properties
	protected static float reactorVolume = 3; // [L]

	protected static float feedFraction = 0.5f; // [dimensionless]

	// duration of the SBR cycle
	protected static float cycleTime = 3; // [h]

	// The SBR sycle
	protected static SbrBulkConcentration.SbrCycle _cycle = new SbrBulkConcentration.SbrCycle(
			cycleTime);

	private static int _nCycles = 400; // write full data concerning cycle

	// every

	// _nCycles

	// Feed composition
	// protected static float oxygenBulkConcentration = 10e-3f; //[gO2/L] DO 100
	protected static float oxygenBulkConcentration = 4e-3f; // [gO2/L] DO 20

	//protected static float ammoniumFeedConcentration = 71e-3f; // [gN/L]
	protected static float ammoniumFeedConcentration = 50e-3f; // [gN/L]

	protected static float nitriteFeedConcentration = 0; // [gN/L]

	protected static float nitrateFeedConcentration = 0; // [gN/L]

	protected static float substrateFeedConcentration = 396e-3f; // [gCOD/L]

	protected static float phosphateBulkConcentration = 20e-3f; // [gP/L]

	private static float oxygenDiffusivity = 8.3e6f; // [um^2/h]

	private static float substrateDiffusivity = 4e6f; // [um2/h]

	private static float phosphateDiffusivity = 4e6f; // [um2/h] ASSUMED!

	private static float ammoniumDiffusivity = 7.083e6f; // [um^2/h]

	private static float nitrateDiffusivity = 6.6667e6f; // [um^2/h]

	private static float nitriteDiffusivity = 6.6667e6f; // [um^2/h]

	// Leave these values
	protected static float precision = 0.01f;

	protected static float ammoniumPrecision = K_XNH_NH4 * precision; // [gN/L]

	protected static float nitritePrecision = K_XNO_NO2 * precision; // [gN/L]

	protected static float nitratePrecision = KPAONO3 * precision; // [gN/L]

	protected static float substratePrecision = KPAOS * precision; // [gCOD/L]

	protected static float phosphatePrecision = KPAOPO4_PP * precision; // [gP/L]

	protected static float maxFractionToDecrease = 0.95f;

	//
	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	protected static float specificMassPolymers = specificMassBiomass * 1e5f; // [gCOD-PHB/L]

	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 1600; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 5.6f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 10 / systemSize;

	protected static float maximumGranuleRadius = 550; // [um]

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side

	// Don't change this

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this

	// detachment rate
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// leave at zero to form round granules

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumberXNH = 20;

	protected static int initialParticleNumberXNO = 20;

	protected static int initialParticleNumberXH = 20;

	protected static int initialParticleNumberXPAO = 20;

	// iteration finish time
	protected static float simulationFinishTime = 17520f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 15.0f; // [h]

	// Computational volume multiplier
	// protected static float nComp = 25.93f * 6.7e5f; //[dimensionless]
	protected static float nComp = 25.93f * 6e5f; // [dimensionless]

	// factor 25.93f makes substrate loading per biomass the same as in real
	// reactor

	// Ammonia and nitrite must be accessible from classes that extend
	// this one so that the controller may be created
	protected SoluteSpecies ammonium;

	protected SoluteSpecies nitrite;

	// /END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// ammonium (NH3)
		ammonium = new SoluteSpecies("ammonium", ammoniumDiffusivity);
		// set up the bulk concentration
		ammonium.setBulkConcentration(new SbrBulkConcentration(
				ammoniumFeedConcentration, feedFraction, ammoniumPrecision,
				maxFractionToDecrease, _cycle));
		// nitrite (NO2)
		nitrite = new SoluteSpecies("nitrite", nitriteDiffusivity);
		// set up the bulk concentration
		nitrite.setBulkConcentration(new SbrBulkConcentration(
				nitriteFeedConcentration, feedFraction, nitritePrecision,
				maxFractionToDecrease, _cycle));
		// nitrite (NO3)
		SoluteSpecies nitrate = new SoluteSpecies("nitrate", nitrateDiffusivity);
		// set up the bulk concentration
		nitrate.setBulkConcentration(new SbrBulkConcentration(
				nitriteFeedConcentration, feedFraction, nitratePrecision,
				maxFractionToDecrease, _cycle));
		// substrate (S)
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		// oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		// set up the bulk concentration for oxygen, uncomment one of these 3:
		// full aeration
		oxygen.setBulkConcentration(createOxygen());
		// Anaerobic feed phase
		// oxygen.setBulkConcentration(new ItermittentBulkConcentration(
		// oxygenBulkConcentration, 0, 2f, 1f));
		// Anaerobic feed phase with controlled shutdown
		// oxygen.setBulkConcentration(new
		// ItermittentBulkConcentrationControlled(
		// oxygenBulkConcentration, 0, 2f, 1f, ammonium, nitrite,
		// nConcentrationsPrecision));
		// set up the bulk concentration
		substrate.setBulkConcentration(new SbrBulkConcentration(
				substrateFeedConcentration, feedFraction, substratePrecision,
				maxFractionToDecrease, _cycle));
		// phosphate (PO4)
		SoluteSpecies phosphate = new SoluteSpecies("phosphate",
				phosphateDiffusivity);
		// set up the bulk concentration
		phosphate.setBulkConcentration(new SbrBulkConcentration(
				phosphateBulkConcentration, feedFraction, phosphatePrecision,
				maxFractionToDecrease, _cycle));
		// 2. Create the particulate species (solids)
		// NH
		Color nhColor = Color.green;
		ParticulateSpecies activeNH = new ParticulateSpecies("activeNH",
				specificMassBiomass, nhColor);
		ParticulateSpecies inertNH = new ParticulateSpecies("inertNH",
				specificMassBiomass, nhColor);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spNH = { activeNH, inertNH };
		float[] fractionalVolumeCompositionNH = { 1.0f, 0 };
		// 3. Create the biomass species
		BiomassSpecies speciesNH = new BiomassSpecies("speciesNH", spNH,
				fractionalVolumeCompositionNH);
		speciesNH.setActiveMass(activeNH);
		speciesNH.setInertMass(inertNH);
		// NHO
		ParticulateSpecies activeNO = new ParticulateSpecies("activeNO",
				specificMassBiomass, Color.yellow);
		ParticulateSpecies inertNO = new ParticulateSpecies("inertNO",
				specificMassBiomass, Color.yellow);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spNO = { activeNO, inertNO };
		float[] fractionalVolumeCompositionNO = { 1.0f, 0 };
		// 3. Create the biomass species
		BiomassSpecies speciesNO = new BiomassSpecies("speciesNO", spNO,
				fractionalVolumeCompositionNO);
		speciesNO.setActiveMass(activeNO);
		speciesNO.setInertMass(inertNO);
		// XH
		ParticulateSpecies activeH = new ParticulateSpecies("activeH",
				specificMassBiomass, Color.blue);
		ParticulateSpecies inertH = new ParticulateSpecies("inertH",
				specificMassBiomass, Color.blue);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spH = { activeH, inertH };
		float[] fractionalVolumeCompositionH = { 1.0f, 0 };
		// 3. Create the biomass species
		BiomassSpecies speciesH = new BiomassSpecies("speciesH", spH,
				fractionalVolumeCompositionH);
		speciesH.setActiveMass(activeH);
		speciesH.setInertMass(inertH);
		// Pao (PAO)
		Color paoColor = Color.red;
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
		ParticulateSpecies[] spPAO = { activePAO, phbPAO, polypPAO,
				glycogenPAO, inertPAO };
		float initialFractionPAO = 0.4f / specificMassBiomass;
		float initialFractionPhb = 0.4f / specificMassPolymers;
		float initialFractionPolyP = 0.1f / specificMassPolymers;
		float initialFractionGlyP = 0.1f / specificMassPolymers;
		float total = initialFractionPAO + initialFractionPolyP
				+ initialFractionGlyP + initialFractionPhb;
		float[] fractionalVolumeCompositionPAO = { initialFractionPAO / total,
				initialFractionPhb / total, initialFractionPolyP / total,
				initialFractionGlyP / total, 0 };
		// 3. Create the biomass species
		BiomassSpecies speciesPAO = new BiomassSpecies("speciesPAO", spPAO,
				fractionalVolumeCompositionPAO);
		speciesPAO.setActiveMass(activePAO);
		speciesPAO.setInertMass(inertPAO);
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// ProcessFactor mS = new Saturation(oxygen, KO);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		// for NH
		ProcessFactor mNHO2 = new Saturation(oxygen, K_XNH_O2);
		ProcessFactor mNHNH4 = new Saturation(ammonium, K_XNH_NH4);
		ProcessFactor iNHNH4 = new Inhibition(ammonium, K_XNH_NH4);
		// //for NO
		ProcessFactor mNOO2 = new Saturation(oxygen, K_XNO_O2);
		ProcessFactor mNONO2 = new Saturation(nitrite, K_XNO_NO2);
		ProcessFactor iNONO2 = new Inhibition(nitrite, K_XNO_NO2);
		// for XH
		ProcessFactor mHO2 = new Saturation(oxygen, K_XH_O2);
		ProcessFactor mHS = new Saturation(substrate, K_XH_S);
		ProcessFactor iHS = new Inhibition(substrate, K_XH_S);
		ProcessFactor iHO2 = new Inhibition(oxygen, K_XH_O2);
		ProcessFactor mHNO3 = new Saturation(nitrate, K_XH_NO3);
		// for PAO
		ProcessFactor mPAOS = new Saturation(substrate, KPAOS);
		ProcessFactor iPAOO2 = new Inhibition(oxygen, KPAOO2);
		ProcessFactor mPAONO3 = new Saturation(nitrate, KPAONO3);
		ProcessFactor mPAONO2 = new Saturation(nitrite, KPAONO2);
		ProcessFactor mf_PAOGLY = new SaturationFromFraction(glycogenPAO,
				activePAO, 0.1f);
		ProcessFactor mf_PAOPP = new SaturationFromFraction(polypPAO,
				activePAO, 0.1f);
		ProcessFactor mPAOO2 = new Saturation(oxygen, KPAOO2);
		// ProcessFactor mPAOPO4 = new Saturation(phosphate, KPAOPO4);
		ProcessFactor iPAOPO4 = new Inhibition(phosphate, KPAOPO4);
		ProcessFactor mPAOPO4_PP = new Saturation(phosphate, KPAOPO4_PP);
		ProcessFactor mf_PAOPHB = new SaturationFromFraction(phbPAO, activePAO,
				KfPHA_P);
		ProcessFactor if_PAOPHB = new InhibitionFromFractionCapacity(phbPAO,
				activePAO, KPHA_P, fPHBmaxPAO);
		ProcessFactor iPAOPHB = new InhibitionFromFraction(phbPAO, activePAO,
				KfPHA_P);
		ProcessFactor iPAOPP = new InhibitionFromFractionCapacity(polypPAO,
				activePAO, KPP_P, fPPmaxPAO);
		ProcessFactor iPAOGly = new InhibitionFromFractionCapacity(glycogenPAO,
				activePAO, 0.01f, fGLYmaxPAO);
		// for all organisms
		// ProcessFactor mNH4 = new Saturation(ammonium, KNH4);
		// 5. Create the reactions
		// NH
		// growth NH
		Reaction growthNH = new Reaction("growthNH", activeNH, uMaxXNH, 2);
		growthNH.addFactor(mNHO2);
		growthNH.addFactor(mNHNH4);
		// NH decay
		Reaction decayNH = new Reaction("decayNH", activeNH, kDecay, 1);
		decayNH.addFactor(iNHNH4);
		// NO
		// growth NO
		Reaction growthNO = new Reaction("growthNO", activeNO, uMaxNO, 2);
		growthNO.addFactor(mNOO2);
		growthNO.addFactor(mNONO2);
		// //growthNO.addFactor(mNH4); //
		// NO decay
		Reaction decayNO = new Reaction("decayNO", activeNO, kDecay, 1);
		decayNO.addFactor(iNONO2);
		// H
		// aerobic growth H
		Reaction aerobicGrowthH = new Reaction("growthH", activeH, uMaxH, 2);
		aerobicGrowthH.addFactor(mHO2);
		aerobicGrowthH.addFactor(mHS);
		// anoxic growth H
		Reaction anoxicGrowthH = new Reaction("anoxicGrowthH", activeH, uMaxH
				* etaNO3, 3);
		anoxicGrowthH.addFactor(iHO2);
		anoxicGrowthH.addFactor(mHNO3);
		anoxicGrowthH.addFactor(mHS);
		// H decay
		Reaction decayH = new Reaction("decayH", activeH, kDecay, 1);
		decayH.addFactor(iHS);
		// PAO
		// substrate uptake
		Reaction sUptakePAO = new Reaction("sUptakePAO", activePAO, qSMaxPAO, 4);
		sUptakePAO.addFactor(mPAOS);
		sUptakePAO.addFactor(if_PAOPHB);
		sUptakePAO.addFactor(mf_PAOGLY);
		sUptakePAO.addFactor(mf_PAOPP);
		// aerobic consumption of PHB
		Reaction aerobicPhbConsumptionPAO = new Reaction(
				"aerobicPhbConsumptionPAO", activePAO, kPhbPAO, 2);
		aerobicPhbConsumptionPAO.addFactor(mf_PAOPHB);
		// aerobicPhbConsumptionPAO.addFactor(mPAOPO4);
		aerobicPhbConsumptionPAO.addFactor(mPAOO2);
		// aerobic storage of PolyP
		Reaction aerobicStoragePolyPPAO = new Reaction(
				"aerobicStoragePolyPPAO", activePAO, kPPPAO, 3);
		aerobicStoragePolyPPAO.addFactor(iPAOPP);
		aerobicStoragePolyPPAO.addFactor(mPAOPO4_PP);
		aerobicStoragePolyPPAO.addFactor(mPAOO2);
		// aerobic storage of glycogen
		Reaction aerobicStorageGlyPAO = new Reaction("aerobicStorageGlyPAO",
				activePAO, kGlyPAO, 3);
		aerobicStorageGlyPAO.addFactor(iPAOGly);
		aerobicStorageGlyPAO.addFactor(mPAOO2);
		aerobicStorageGlyPAO.addFactor(mf_PAOPHB);
		// anoxic (NO3) consumption of PHB
		Reaction anoxicNO3PhbConsumptionPAO = new Reaction(
				"anoxicNO3PhbConsumptionPAO", activePAO, kPhbPAO * etaNO3, 3);
		anoxicNO3PhbConsumptionPAO.addFactor(mf_PAOPHB);
		// anoxicNO3PhbConsumptionPAO.addFactor(mPAOPO4);
		anoxicNO3PhbConsumptionPAO.addFactor(iPAOO2);
		anoxicNO3PhbConsumptionPAO.addFactor(mPAONO3);
		// anoxic (NO3) storage of PolyP
		Reaction anoxicNO3StoragePolyPPAO = new Reaction(
				"anoxicNO3StoragePolyPPAO", activePAO, kPPPAO * etaNO3, 4);
		anoxicNO3StoragePolyPPAO.addFactor(iPAOPP);
		anoxicNO3StoragePolyPPAO.addFactor(mPAOPO4_PP);
		anoxicNO3StoragePolyPPAO.addFactor(iPAOO2);
		anoxicNO3StoragePolyPPAO.addFactor(mPAONO3);
		// anoxic (NO3) storage of glycogen
		Reaction anoxicNO3StorageGlyPAO = new Reaction(
				"anoxicNO3StorageGlyPAO", activePAO, kGlyPAO * etaNO3, 4);
		anoxicNO3StorageGlyPAO.addFactor(iPAOGly);
		anoxicNO3StorageGlyPAO.addFactor(mf_PAOPHB);
		anoxicNO3StorageGlyPAO.addFactor(iPAOO2);
		anoxicNO3StorageGlyPAO.addFactor(mPAONO3);
		// anoxic (NO2) consumption of PHB
		Reaction anoxicNO2PhbConsumptionPAO = new Reaction(
				"anoxicNO2PhbConsumptionPAO", activePAO, kPhbPAO * etaNO2, 3);
		anoxicNO2PhbConsumptionPAO.addFactor(mf_PAOPHB);
		// anoxicNO2PhbConsumptionPAO.addFactor(mPAOPO4);
		anoxicNO2PhbConsumptionPAO.addFactor(iPAOO2);
		anoxicNO2PhbConsumptionPAO.addFactor(mPAONO2);
		// anoxic (NO2) storage of PolyP
		Reaction anoxicNO2StoragePolyPPAO = new Reaction(
				"anoxicNO2StoragePolyPPAO", activePAO, kPPPAO * etaNO2, 4);
		anoxicNO2StoragePolyPPAO.addFactor(iPAOPP);
		anoxicNO2StoragePolyPPAO.addFactor(mPAOPO4_PP);
		anoxicNO2StoragePolyPPAO.addFactor(iPAOO2);
		anoxicNO2StoragePolyPPAO.addFactor(mPAONO2);
		// anoxic (NO2) storage of glycogen
		Reaction anoxicNO2StorageGlyPAO = new Reaction(
				"anoxicNO2StorageGlyPAO", activePAO, kGlyPAO * etaNO2, 4);
		anoxicNO2StorageGlyPAO.addFactor(iPAOGly);
		anoxicNO2StorageGlyPAO.addFactor(mf_PAOPHB);
		anoxicNO2StorageGlyPAO.addFactor(iPAOO2);
		anoxicNO2StorageGlyPAO.addFactor(mPAONO2);
		// decay
		Reaction decayPAO = new Reaction("decayPAO", activePAO, kDecay, 2);
		decayPAO.addFactor(iPAOPHB);
		decayPAO.addFactor(iPAOPO4);
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass NH
		NetReaction rsNHactive = new NetReaction(2);
		rsNHactive.addReaction(growthNH, 1);
		rsNHactive.addReaction(decayNH, -1);
		activeNH.setProcesses(rsNHactive);
		// inert NH
		NetReaction rsInertNH = new NetReaction(1);
		rsInertNH.addReaction(decayNH, 1);
		inertNH.setProcesses(rsInertNH);
		// active NO
		NetReaction rsActiveNO = new NetReaction(2);
		rsActiveNO.addReaction(growthNO, 1);
		rsActiveNO.addReaction(decayNO, -1);
		activeNO.setProcesses(rsActiveNO);
		// inert NO
		NetReaction rsInertNO = new NetReaction(1);
		rsInertNO.addReaction(decayNO, fIDecay);
		inertNO.setProcesses(rsInertNO);
		// active H
		NetReaction rsActiveH = new NetReaction(3);
		rsActiveH.addReaction(aerobicGrowthH, 1);
		rsActiveH.addReaction(anoxicGrowthH, 1);
		rsActiveH.addReaction(decayH, -1);
		activeH.setProcesses(rsActiveH);
		// inert H
		NetReaction rsInertH = new NetReaction(1);
		rsInertH.addReaction(decayH, fIDecay);
		inertH.setProcesses(rsInertH);
		// active mass PAO
		NetReaction rsPaoActive = new NetReaction(10);
		rsPaoActive.addReaction(aerobicPhbConsumptionPAO, 1 / YPhbAerobicPAO);
		rsPaoActive.addReaction(anoxicNO3PhbConsumptionPAO,
				1 / YPhbAnoxicNO3PAO);
		rsPaoActive.addReaction(anoxicNO2PhbConsumptionPAO,
				1 / YPhbAnoxicNO2PAO);
		rsPaoActive.addReaction(aerobicStorageGlyPAO, -1 / YGlyAerobicPAO);
		rsPaoActive.addReaction(anoxicNO3StorageGlyPAO, -1 / YGlyAnoxicNO3PAO);
		rsPaoActive.addReaction(anoxicNO2StorageGlyPAO, -1 / YGlyAnoxicNO2PAO);
		rsPaoActive.addReaction(aerobicStoragePolyPPAO, -1 / YPpAerobicPAO);
		rsPaoActive.addReaction(anoxicNO3StoragePolyPPAO, -1 / YPpAnoxicNO3PAO);
		rsPaoActive.addReaction(anoxicNO2StoragePolyPPAO, -1 / YPpAnoxicNO2PAO);
		rsPaoActive.addReaction(decayPAO, -1);
		activePAO.setProcesses(rsPaoActive);
		// PHB in PAO
		NetReaction rsPhbPAO = new NetReaction(4);
		rsPhbPAO.addReaction(sUptakePAO, Y_PHB_S_of_PAO);
		rsPhbPAO.addReaction(aerobicPhbConsumptionPAO, -1);
		rsPhbPAO.addReaction(anoxicNO3PhbConsumptionPAO, -1);
		rsPhbPAO.addReaction(anoxicNO2PhbConsumptionPAO, -1);
		phbPAO.setProcesses(rsPhbPAO);
		// PP in PAO
		NetReaction rsPpPAO = new NetReaction(4);
		rsPpPAO.addReaction(sUptakePAO, -Y_PO4_S);
		rsPpPAO.addReaction(aerobicStoragePolyPPAO, 1);
		rsPpPAO.addReaction(anoxicNO3StoragePolyPPAO, 1);
		rsPpPAO.addReaction(anoxicNO2StoragePolyPPAO, 1);
		polypPAO.setProcesses(rsPpPAO);
		// Glycogen in PAO
		NetReaction rsGlyPAO = new NetReaction(4);
		rsGlyPAO.addReaction(sUptakePAO, 1 - Y_PHB_S_of_PAO);
		rsGlyPAO.addReaction(aerobicStorageGlyPAO, 1);
		rsGlyPAO.addReaction(anoxicNO3StorageGlyPAO, 1);
		rsGlyPAO.addReaction(anoxicNO2StorageGlyPAO, 1);
		glycogenPAO.setProcesses(rsGlyPAO);
		// inerts in PAO
		NetReaction rsInertPAO = new NetReaction(1);
		rsInertPAO.addReaction(decayPAO, 1);
		inertPAO.setProcesses(rsInertPAO);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(6);
		rsOxygen.addReaction(growthNH, 1 - 3.43f / Y_XNH_NH4);
		rsOxygen.addReaction(growthNO, 1 - 1.14f / Y_XNO_NO2);
		rsOxygen.addReaction(aerobicGrowthH, 1 - 1f / Y_XH_S);
		rsOxygen.addReaction(aerobicPhbConsumptionPAO, 1 / YPhbAerobicPAO - 1);
		rsOxygen.addReaction(aerobicStorageGlyPAO, 1 - 1 / YGlyAerobicPAO);
		rsOxygen.addReaction(aerobicStoragePolyPPAO, -(1 / YPpAerobicPAO));
		oxygen.setProcesses(rsOxygen);
		// ammonium
		NetReaction rsAmmonium = new NetReaction(1);
		rsAmmonium.addReaction(growthNH, -1 / Y_XNH_NH4);
		ammonium.setProcesses(rsAmmonium);
		// nitrite (NO2)
		NetReaction rsNitrite = new NetReaction(5);
		rsNitrite.addReaction(growthNH, 1 / Y_XNH_NH4);
		rsNitrite.addReaction(growthNO, -1 / Y_XNO_NO2);
		rsNitrite.addReaction(anoxicNO2PhbConsumptionPAO,
				(1 / YPhbAnoxicNO2PAO - 1) / 1.71f);
		rsNitrite.addReaction(anoxicNO2StorageGlyPAO,
				(1 - 1 / YGlyAnoxicNO2PAO) / 1.71f);
		rsNitrite.addReaction(anoxicNO2StoragePolyPPAO,
				-(1 / 1.71f / YPpAnoxicNO3PAO));
		nitrite.setProcesses(rsNitrite);
		// nitrate (NO3)
		NetReaction rsNitrate = new NetReaction(5);
		rsNitrate.addReaction(growthNO, 1 / Y_XNO_NO2);
		rsNitrate.addReaction(anoxicGrowthH, (1 - 1f / Y_XH_S) / 2.86f);
		rsNitrate.addReaction(anoxicNO3PhbConsumptionPAO,
				(1 / YPhbAnoxicNO3PAO - 1) / 2.86f);
		rsNitrate.addReaction(anoxicNO3StoragePolyPPAO,
				-(1 / 2.86f / YPpAnoxicNO3PAO));
		rsNitrate.addReaction(anoxicNO3StorageGlyPAO,
				(1 - 1 / YGlyAnoxicNO3PAO) / 2.86f);
		nitrate.setProcesses(rsNitrate);
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(3);
		rsSubstrate.addReaction(aerobicGrowthH, -1 / Y_XH_S);
		rsSubstrate.addReaction(anoxicGrowthH, -1 / Y_XH_S);
		rsSubstrate.addReaction(sUptakePAO, -1);
		substrate.setProcesses(rsSubstrate);
		// phosphate (PO4)
		NetReaction rsPhosphate = new NetReaction(3);
		rsPhosphate.addReaction(sUptakePAO, Y_PO4_S);
		rsPhosphate.addReaction(aerobicStoragePolyPPAO, -1);
		rsPhosphate.addReaction(anoxicNO3StoragePolyPPAO, -1);
		phosphate.setProcesses(rsPhosphate);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesNH);
		addBiomassSpecies(speciesNO);
		addBiomassSpecies(speciesH);
		addBiomassSpecies(speciesPAO);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(ammonium);
		addSoluteSpecies(nitrite);
		addSoluteSpecies(nitrate);
		addSoluteSpecies(substrate);
		addSoluteSpecies(phosphate);
	}

	/**
	 * @return the bulk concentration of oxygen
	 */
	public abstract BulkConcentration createOxygen();

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumberXNH, initialParticleNumberXNO,
				initialParticleNumberXH, initialParticleNumberXPAO };
		inoculateRandomly(nCells);
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
			Model.model().setVerticalCutoffSize(maximumGranuleRadius);
		} catch (InvalidValueException e) {
			System.out.println(e);
			System.exit(-1);
		}
	}

	/**
	 * Set the reactors dimensions
	 */
	private static void setTheReactorDimensions(ApplicationComponent app) {
		float rVIM = reactorVolume * 1e15f; // reactor volume in cubic
		// micrometer
		float carrierArea = nComp
				* Model.model().getComputationalVolumeCarrierArea();
		app.setReactorParameters(cycleTime, carrierArea, rVIM);
	}

	/**
	 * @param app
	 */
	protected static void setSystemParametersAndInitializeSystemSpace(
			ApplicationComponent app) {
		// create the space
		app.setSystemSpaceParameters(geometry, systemSize,
				relativeMaximumRadius, relativeMinimumRadius,
				relativeBoundaryLayer, gridSide, kShov);
		// initialize
		app.initializeSystemSpace();
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public void run(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		boolean runWithGraphics = false;
		if (args.length != 2) {
			throw new RuntimeException("Program arguments missing: "
					+ "2 program arguments should be supplied"
					+ " (1 - output directory,"
					+ " 2 - flag for running with graphics)");
		}
		// parse the input arguments
		outputDirectory = args[0];
		int arg1 = Integer.parseInt(args[1]);
		switch (arg1) {
		case 0:
			runWithGraphics = false;
			break;
		case 1:
			runWithGraphics = true;
			break;
		default:
			throw new RuntimeException("second program" + " argument must be 0"
					+ " (for running with no graphics) "
					+ "or 1 (for running with graphics)");
		}
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = this;
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VariableSeries biovolume = new BiovolumeSeries();
		VariableSeries[] runLengthSeries = { new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries() };
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// bulk concentrations
			app = new BulkConcentrationVizualizer(app);
			// finally, the controller must be the last decorator to add
			app = new VizualModelControler(app);
		}
		try {
			// create the space
			setSystemParametersAndInitializeSystemSpace(app);
			// initialize
			app.intializeStateWriters(outputDirectory);
			// Pov witer is added twice
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
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(runLengthSeries[0]);
			spw.addSeries(runLengthSeries[1]);
			spw.addSeries(runLengthSeries[2]);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// add the time constraints writer
			app.addStateWriter(new TimeConstraintsWriter());
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
		} catch (ModelException e) {
			System.out.println(e);
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
			// start the iteration
			app.writeState(); // write iteration 0
			app.startIterating();
		} catch (Exception e1) {
			try {
				app.forceWriteState();
			} catch (IOException e2) {
				System.err.println("Error serializing state:");
				System.err.println("");
				e2.printStackTrace();
			}
			System.err.println("");
			System.err.println("Program failed due to :");
			System.err.println("");
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}