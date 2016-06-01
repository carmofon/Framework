package nl.tudelft.bt.model.work.will;

import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.work.benzoate.MotileBiomassSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.reaction.*;

public class ChemotacticBiomassSpecies extends MotileBiomassSpecies {
	
	private SoluteSpecies _attractant;
	
	private float _Kspeed;

	private float _Kspread;
	
	private float _diffusivity;
	
	protected ProcessFactor _motilityFactor;
	
	protected float _minimumMotility;
	
	/**
	 * 
	 * @param name
	 * 				the name of the chemotactic biomass species
	 * @param species
	 * 				array of constituent particulate species
	 * @param fractionalCompositionInVolume
	 * 				array of fractional compositions for particulates
	 * @param attractant
	 * 				solute species to which this biomass is attracted
	 * @param motilityDiffusivity
	 * 				the max motility of this species
	 * @param Kspeed
	 * 				half-saturation constant for motility of this species
	 * @param Kspread
	 * 				half-saturation constant for chemotactic directedness of this species
	 * @throws NonMatchingNumberException
	 */
	public ChemotacticBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, SoluteSpecies attractant, float motilityDiffusivity,
			float Kspeed, float Kspread, float minMotility)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume, motilityDiffusivity);
		_attractant = attractant;
		_Kspeed = Kspeed;
		_Kspread = Kspread;
		_diffusivity = motilityDiffusivity;
		_minimumMotility = minMotility;
		_motilityFactor = new SaturationWithFloor(_attractant, _Kspeed, minMotility);
//		_motilityFactor = new Hill(_attractant, (float) Math.pow(_Kspeed, 2), 3);
//		_motilityFactor = new Step(_attractant, _Kspeed);
	}
	
	public BiomassParticle createBiomassParticle() {
		return new ChemotacticBiomassParticle(this);
	}
	
	public SoluteSpecies getChemoattractant () {
		return _attractant;
	}
	
	public float getKspeed() {
		return _Kspeed;
	}
	
	public float getKspread() {
		return _Kspread;
	}
	
	public float getDiffusivity() {
		return _diffusivity;
	}

}
