/**
 * 
 */
package nl.tudelft.bt.model.work.will.cell;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.work.will.MutatorBiomassSpecies;
import nl.tudelft.bt.model.work.will.PolarizableBiomassParticle;

/**
 * @author wkc
 *
 */
public class ZonalVicsekBiomassSpecies extends MutatorBiomassSpecies {

	// stop simulation on periodic boundary crossing?
	protected boolean _stopOnPeriodicBoundaryCrossing = false;

	// chemoattractant species
	protected SoluteSpecies _chemoattractant;
	
	// gradient sensing parameter
	protected float _wg;
	
	// neighbor polarization sensing parameter
	protected float _ws;
	
	// random polarization parameter
	protected float _wr;
	
	// motility diffusivity
	protected float _motility;
	
	// relative adhesive radius
	protected float _ra;
	
	// sensitivity to chemoattractant concentration
	protected float _kc;
	
	// sensitivity to chemoattractant gradient
	protected float _kg;
	
	// chemoattractant threshold
//	protected float _cth;
	
	// Hill coefficient for tuning response
	protected float _co = 2;
	
	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param wildTypeTrait
	 * @param mutantTrait
	 * @param mutationRate
	 * @throws NonMatchingNumberException
 	 * @param chemoattractant	the chemoattractant solute species
	 * @param motility	the motility diffusivity
	 * @param RA	the relative adhesive radius
	 * @param wg	the weight of gradient sensing signals
	 * @param ws	the weight of cell-cell coordination signals
	 * @param wr	the weight of random polarization
	 * @param kc	the response constant for chemoattractant sensitivity
	 * @param kg	the response constant for gradient sensitivity
	 */
	public ZonalVicsekBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume,
			ParticulateSpecies wildTypeTrait, ParticulateSpecies mutantTrait,
			float mutationRate, SoluteSpecies chemoattractant,
			float motility, float RA, float wg, float ws, float wr,
			float kc, float kg) throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume, wildTypeTrait,
				mutantTrait, mutationRate);
		_chemoattractant = chemoattractant;
		_motility = motility;
		_wg = wg;
		_ws = ws;
		_wr = wr;
		_ra = RA;
		_kc = kc;
		_kg = kg;
	}
	
	@Override
	public BiomassParticle createBiomassParticle() {
		return new ZonalVicsekBiomassParticle(this);
	}
	
	/**
	 * Set simulation to stop when a cell of this species crosses periodic
	 * boundary
	 */
	public void setStopOnPeriodicBoundaryCrossing() {
		_stopOnPeriodicBoundaryCrossing = true;
	}

}